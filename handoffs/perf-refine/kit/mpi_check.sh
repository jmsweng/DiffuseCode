#!/bin/bash
# Handoff 09: do refinements under MPI give the same result as serial ones?
# Builds an MPI-enabled discus_suite from the tip of origin/perf/refine, and a variant with the old
# step-size bug (kit/mpi_delta_unfix.patch). Runs four 4-parameter fits serially (bins/c7 of handoff 00)
# and under mpirun with 5 ranks (4 workers, so every worker is used in every cycle), and compares them.
#   bash mpi_check.sh
# Environment overrides (to reuse prebuilt binaries): MPI_TIP, MPI_UNFIXED, SERIAL_BIN, MPIRUN_OPTS
set -uo pipefail
ROOT=~/refine_validation
KIT=$ROOT/kit
M=$ROOT/mpi
SRC=$ROOT/src
SERIAL_BIN=${SERIAL_BIN:-$ROOT/bins/c7}
MPI_TIP=${MPI_TIP:-$M/mpi_tip}
MPI_UNFIXED=${MPI_UNFIXED:-$M/mpi_tip_unfixed}
mkdir -p $M

fail() { echo "RESULT: FAIL ($1)"; exit 1; }

command -v mpif90 >/dev/null && command -v mpirun >/dev/null && command -v script >/dev/null \
  || fail "OpenMPI or 'script' missing; needs: sudo apt install -y libopenmpi-dev openmpi-bin bsdutils"
[ -x "$SERIAL_BIN" ] || fail "$SERIAL_BIN missing (run handoff 00 first)"

# ---- 1. build the two MPI binaries (skipped when they exist) -------------------------------
if [ ! -x "$MPI_TIP" ] || [ ! -x "$MPI_UNFIXED" ]; then
  git -C $SRC fetch -q origin || fail "git fetch"
  TIP=$(git -C $SRC rev-parse origin/perf/refine) || fail "origin/perf/refine not found"
  [ -d $M/wt ] || git -C $SRC worktree add --detach $M/wt $TIP >/dev/null 2>&1 || fail "git worktree"
  git -C $M/wt checkout -q -- . && git -C $M/wt checkout -q --detach $TIP || fail "checkout $TIP"
  mkdir -p $M/build && cd $M/build
  if [ ! -f CMakeCache.txt ]; then
    FC=$(which gfortran) PGPLOT_DIR=${PGPLOT_DIR:-/usr/lib/pgplot5} \
      cmake $M/wt -DDIFFUSE_GRAPHICS=OFF -DCMAKE_BUILD_TYPE=Release > cmake.log 2>&1 || fail "cmake, see $M/build/cmake.log"
  fi
  grep -q "building SUITE with MPI" cmake.log && grep -q "building DIFFEV with MPI" cmake.log \
    || fail "cmake did not find MPI; see $M/build/cmake.log"
  echo "building MPI binary from $TIP ..."
  make -j$(nproc) > make_tip.log 2>&1 || fail "make, see $M/build/make_tip.log"
  cp suite/prog/discus_suite $M/mpi_tip
  git -C $M/wt apply $KIT/mpi_delta_unfix.patch || fail "kit/mpi_delta_unfix.patch does not apply"
  echo "building the variant with the old step-size bug ..."
  make -j$(nproc) > make_unfixed.log 2>&1 || fail "make (unfixed), see $M/build/make_unfixed.log"
  cp suite/prog/discus_suite $M/mpi_tip_unfixed
  git -C $M/wt checkout -q -- .   # leave the scratch source clean
fi
[ -x "$MPI_TIP" ] && [ -x "$MPI_UNFIXED" ] || fail "MPI binaries missing"

# ---- 2. test cases ---------------------------------------------------------------------------
# The work macro is self-contained (MPI workers have no data loaded) and has no counters.
T=$M/test_$(date +%H%M%S)
mkdir -p $T
cat > $T/work.mac <<'EOF'
branch kuplot
func P_amp*exp(0.0-((r[0]-P_cen)/P_wid)**2) + P_bg, -5, 5, 0.1
exit
finished
EOF
CASES="det far p5 p2far"   # easy start; far start with rejected steps; 5-point derivatives; points:2 far start
for c in $CASES; do
  for v in serial mpi; do
    mkdir -p $T/$c/$v
    cp $KIT/refine/gauss.data $T/work.mac $T/$c/$v/
    { echo 'set prompt,off'
      echo 'seed 4711, 271828, 314159, 161803, 141421, 173205, 223606, 99991'
      echo refine; echo 'data xy, gauss.data'; cat $KIT/refine/p_$c.txt
      echo 'set log, on'; echo 'set cycle, 25'; echo 'run work.mac'; echo exit; echo exit; } > $T/$c/$v/main.mac
  done
  (cd $T/$c/serial && timeout 300 $SERIAL_BIN -macro main.mac < /dev/null > out.txt 2>&1)
  # 'script' gives the program a terminal. Without one it cannot find its parent process
  # (it runs 'ps j'), silently decides MPI is off and runs serially.
  (cd $T/$c/mpi && script -qec "timeout 400 mpirun ${MPIRUN_OPTS:-} --oversubscribe -n 5 $MPI_TIP -macro main.mac < /dev/null" /dev/null > out.txt 2>&1)
done
# the old bug, 5-point case only
mkdir -p $T/p5/unfixed && cp $T/p5/mpi/main.mac $T/p5/mpi/gauss.data $T/work.mac $T/p5/unfixed/
(cd $T/p5/unfixed && script -qec "timeout 400 mpirun ${MPIRUN_OPTS:-} --oversubscribe -n 5 $MPI_UNFIXED -macro main.mac < /dev/null" /dev/null > out.txt 2>&1)

# ---- 3. evaluate -----------------------------------------------------------------------------
OUT=$(python3 - "$T" "$CASES" <<'EOF'
import re, sys, glob
T, cases = sys.argv[1], sys.argv[2].split()
def parse(path):
    t = open(path, errors='replace').read()
    t = re.sub(r'\x1b\[[0-9;]*m', '', t)
    pars = {m.group(1): (float(m.group(2)), float(m.group(3)))
            for m in re.finditer(r'^\s+(P_\w+)\s+:\s+(\S+)\s+\+-\s+(\S+)', t, re.M)}
    cyc = [l.split() for l in t.splitlines()
           if re.match(r'^\s*\d+\s+\S+\s+\S+\s+\d+\s+\S+\s+\S+\s+\S+\s+\S+\s*$', l)]
    errs = sorted({l.strip() for l in t.splitlines() if re.search(r'\*\*\*[A-Z/]+\*\*\* ', l)
                   and 'Version' not in l and '*****' not in l})
    return pars, cyc, errs
ok = True
for c in cases:
    s, m = parse(f'{T}/{c}/serial/out.txt'), parse(f'{T}/{c}/mpi/out.txt')
    logs = len(glob.glob(f'{T}/{c}/mpi/DISCUS_SUITE_DERIVATIVES/LOGFILE.*'))
    same_pars = bool(s[0]) and set(s[0]) == set(m[0]) and all(
        abs(s[0][k][0] - m[0][k][0]) <= 1e-4 * max(1.0, abs(s[0][k][0])) for k in s[0])
    same_cyc = len(s[1]) == len(m[1]) and len(s[1]) > 0
    good = same_pars and same_cyc and not m[2] and logs > 0
    ok &= good
    print(f'{"PASS" if good else "FAIL"}  {c:6s} serial: {len(s[1])} cycles, MPI: {len(m[1])} cycles, parameters {"equal" if same_pars else "DIFFER"}, '
          f'MPI errors {len(m[2])}, worker logs {logs}')
    if not good:
        for e in m[2][:3]:
            print('        MPI error:', e)
        if logs == 0:
            print('        no worker logs: MPI was not active, the run was serial')
# does the case detect the old bug?
s, u = parse(f'{T}/p5/serial/out.txt'), parse(f'{T}/p5/unfixed/out.txt')
c0s = float(s[1][0][1]) if s[1] else None
c0u = float(u[1][0][1]) if u[1] else None
detect = c0s is not None and c0u is not None and abs(c0u - c0s) > 0.01 * abs(c0s)
print(f'{"PASS" if detect else "INFO"}  p5     with the old step-size bug, chi2/(N-P) in cycle 0 is {c0u} (serial {c0s}): '
      f'{"the old bug is detected" if detect else "the old bug is NOT detected by this case"}')
print('RESULT: PASS' if ok and detect else ('RESULT: FAIL' if not ok else 'RESULT: INCONCLUSIVE (the case does not detect the old bug)'))
sys.exit(0 if ok and detect else (1 if not ok else 2))
EOF
)
rc=$?
echo "$OUT" | grep -v '^RESULT:'
if [ $rc -eq 1 ]; then
  for c in $CASES; do
    echo "=== end of output of $c / mpi"
    sed 's/\x1b\[[0-9;]*m//g' $T/$c/mpi/out.txt | grep -v '^\s*$' | tail -n 12
  done
fi
echo "(raw outputs: $T)"
echo "$OUT" | grep '^RESULT:'
exit $rc

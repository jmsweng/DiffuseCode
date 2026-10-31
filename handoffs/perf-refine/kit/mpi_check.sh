#!/bin/bash
# Handoff 09: does the MPI derivative fix (c1) change a real MPI refinement?
# Builds two MPI-enabled binaries (mpi_c0k = upstream refine_run.f90, mpi_c1k = with the fix; both with
# the KUPLOT exit patch mpi_kuplot_exit.patch, without which no MPI fit with a KUPLOT curve can finish),
# runs a 2-parameter, 5-point-derivative fit serially (reference) and under mpirun, and compares the first cycle.
#   bash mpi_check.sh            full run (build if needed)
# Environment overrides (to reuse prebuilt binaries): MPI_C0, MPI_C1, SERIAL_BIN
set -uo pipefail
ROOT=~/refine_validation
KIT=$ROOT/kit
M=$ROOT/mpi
SRC=$ROOT/src
mapfile -t SHAS < <(grep -v '^#' $KIT/COMMITS | awk '{print $1}')
C0=${SHAS[0]}; C1=${SHAS[1]}
MPI_C0=${MPI_C0:-$M/mpi_c0k}
MPI_C1=${MPI_C1:-$M/mpi_c1k}
SERIAL_BIN=${SERIAL_BIN:-$ROOT/bins/c1}
mkdir -p $M

fail() { echo "RESULT: FAIL ($1)"; exit 1; }
inconclusive() { echo "RESULT: INCONCLUSIVE ($1)"; exit 2; }

command -v mpif90 >/dev/null && command -v mpirun >/dev/null && command -v script >/dev/null \
  || fail "OpenMPI or 'script' missing; needs: sudo apt install -y libopenmpi-dev openmpi-bin bsdutils"

# ---- 1. build the two MPI binaries (skipped when they exist) -------------------------------
if [ ! -x "$MPI_C1" ] || [ ! -x "$MPI_C0" ]; then
  [ -d $M/wt ] || git -C $SRC worktree add --detach $M/wt $C1 >/dev/null 2>&1 || fail "git worktree"
  git -C $M/wt checkout -q --detach $C1
  git -C $M/wt checkout -q -- kuplot/prog/kuplot_exit.f90 refine/prog/refine_run.f90
  git -C $M/wt apply $KIT/mpi_kuplot_exit.patch || fail "kit/mpi_kuplot_exit.patch does not apply"
  mkdir -p $M/build && cd $M/build
  if [ ! -f CMakeCache.txt ]; then
    FC=$(which gfortran) PGPLOT_DIR=${PGPLOT_DIR:-/usr/lib/pgplot5} \
      cmake $M/wt -DDIFFUSE_GRAPHICS=OFF -DCMAKE_BUILD_TYPE=Release > cmake.log 2>&1 || fail "cmake, see $M/build/cmake.log"
  fi
  grep -q "building SUITE with MPI" cmake.log && grep -q "building DIFFEV with MPI" cmake.log \
    || fail "cmake did not find MPI; see $M/build/cmake.log"
  echo "building MPI c1 (fixed) ..."
  make -j$(nproc) > make_c1.log 2>&1 || fail "make c1, see $M/build/make_c1.log"
  cp suite/prog/discus_suite $M/mpi_c1k
  # c0 differs from c1 only in refine/prog/refine_run.f90
  git -C $M/wt checkout $C0 -- refine/prog/refine_run.f90
  echo "building MPI c0 (unfixed) ..."
  make -j$(nproc) > make_c0.log 2>&1 || fail "make c0, see $M/build/make_c0.log"
  cp suite/prog/discus_suite $M/mpi_c0k
  git -C $M/wt checkout $C1 -- refine/prog/refine_run.f90
  git -C $M/wt checkout -q -- kuplot/prog/kuplot_exit.f90   # leave the scratch source clean
fi
[ -x "$MPI_C0" ] && [ -x "$MPI_C1" ] && [ -x "$SERIAL_BIN" ] || fail "binaries missing (run handoff 00 first)"

# ---- 2. test case: 2 parameters, 5-point derivatives => 8 derivative runs per cycle ---------
# The work macro is self-contained (the MPI workers have no data loaded) and has no counters.
T=$M/test_$(date +%H%M%S)
mkdir -p $T
mk() {  # dir
  mkdir -p $T/$1
  cp $KIT/refine/gauss.data $T/$1/
  cat > $T/$1/work_mpi.mac <<'EOF'
branch kuplot
func P_amp*exp(0.0-((r[0]-0.7)/1.3)**2) + P_bg, -5, 5, 0.1
exit
finished
EOF
  cat > $T/$1/main.mac <<'EOF'
set prompt,off
seed 4711, 271828, 314159, 161803, 141421, 173205, 223606, 99991
refine
data xy, gauss.data
newparam P_amp, value:6.0, points:5
newparam P_bg,  value:0.2, points:5
set log, on
set cycle, 3
run work_mpi.mac
exit
exit
EOF
}
for d in serial mpi_c0 mpi_c1; do mk $d; done

NRANK=9   # 1 master + 8 workers; one derivative run per worker (see notes in 09_MPI.md)
(cd $T/serial && timeout 300 $SERIAL_BIN -macro main.mac < /dev/null > out.txt 2>&1)
for v in mpi_c0 mpi_c1; do
  bin=$MPI_C0; [ $v = mpi_c1 ] && bin=$MPI_C1
  # 'script' gives the program a terminal. Without one it cannot find its parent process
  # (it runs 'ps j'), silently decides MPI is off and runs serially.
  (cd $T/$v && script -qec "timeout 400 mpirun ${MPIRUN_OPTS:-} --oversubscribe -n $NRANK $bin -macro main.mac < /dev/null" /dev/null > out.txt 2>&1)
done

# ---- 3. evaluate -----------------------------------------------------------------------------
OUT=$(python3 - "$T" <<'EOF'
import re, sys, os, glob
T = sys.argv[1]
def parse(v):
    t = open(f'{T}/{v}/out.txt', errors='replace').read()
    t = re.sub(r'\x1b\[[0-9;]*m', '', t)
    pars = {m.group(1): float(m.group(2)) for m in re.finditer(r'^\s+(P_\w+)\s+:\s+(\S+)\s+\+-', t, re.M)}
    cyc = [l for l in t.splitlines() if re.match(r'^\s*\d+\s+\S+\s+\S+\s+\d+\s+\S+\s+\S+\s+\S+\s+\S+\s*$', l)]
    errs = sorted({l.strip() for l in t.splitlines() if '***' in l and 'Version' not in l and '*****' not in l})
    return pars, cyc, errs
res = {v: parse(v) for v in ('serial', 'mpi_c0', 'mpi_c1')}
def cycle0(c):   # chi^2/(N-P) of the first cycle
    return float(c[0].split()[1]) if c else None
for v, (p, c, e) in res.items():
    print(f'{v:8s} cycles={len(c):2d} chi2/(N-P) of cycle 0 = {cycle0(c)}   final pars={ {k: round(x, 5) for k, x in p.items()} }')
    for x in e[:3]:
        print('         error:', x)
# MPI really on? the workers write LOGFILE.nnnn into DISCUS_SUITE_DERIVATIVES (kept by 'set log, on')
for v in ('mpi_c0', 'mpi_c1'):
    n = len(glob.glob(f'{T}/{v}/DISCUS_SUITE_DERIVATIVES/LOGFILE.*'))
    print(f'{v:8s} worker log files: {n}')
s0, a, b = cycle0(res['serial'][1]), cycle0(res['mpi_c1'][1]), cycle0(res['mpi_c0'][1])
logs1 = len(glob.glob(f'{T}/mpi_c1/DISCUS_SUITE_DERIVATIVES/LOGFILE.*'))
if s0 is None:
    print('RESULT: FAIL (serial reference fit did not run)'); sys.exit(1)
if logs1 == 0:
    print('RESULT: FAIL (MPI was not active: no worker logs; the run was serial)'); sys.exit(1)
if a is None:
    print('RESULT: INCONCLUSIVE (the fixed MPI fit has no cycle 0; see the errors above)'); sys.exit(2)
if abs(a - s0) > 1e-3 * abs(s0):
    print('RESULT: FAIL (cycle 0 of the fixed MPI fit differs from the serial fit: %s against %s)' % (a, s0)); sys.exit(1)
if b is not None and abs(b - s0) <= 1e-2 * abs(s0):
    print('RESULT: INCONCLUSIVE (fixed MPI equals serial, but the unfixed build also does: the case cannot detect the bug)'); sys.exit(2)
print('PASS  cycle 0: serial %s, MPI with the fix %s, MPI without the fix %s' % (s0, a, b))
print('NOTE  later cycles can show "Error opening file" (a worker writing its result file a second time): known, not checked here')
print('RESULT: PASS')
EOF
)
rc=$?
echo "$OUT" | grep -v '^RESULT:'
if [ $rc -ne 0 ]; then
  for v in serial mpi_c0 mpi_c1; do
    echo "=== end of output of $v"
    sed 's/\x1b\[[0-9;]*m//g' $T/$v/out.txt | grep -v '^\s*$' | tail -n 20
  done
fi
echo "(raw outputs: $T)"
echo "$OUT" | grep '^RESULT:'
exit $rc

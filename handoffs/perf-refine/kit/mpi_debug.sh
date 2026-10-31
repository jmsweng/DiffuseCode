#!/bin/bash
# Handoff 10: show what the MPI master reads after the workers return (information gathering).
# Needs handoff 09 to have built ~/refine_validation/mpi (worktree 'wt' at c1 and build dir 'build').
# The debug build is made from a scratch copy of the source in ~/refine_validation/mpi; the repository is not touched.
set -uo pipefail
ROOT=~/refine_validation
KIT=$ROOT/kit
M=$ROOT/mpi
fail() { echo "RESULT: FAIL ($1)"; exit 1; }
[ -d $M/wt ] && [ -f $M/build/CMakeCache.txt ] || fail "run handoff 09 first (it creates $M/wt and $M/build)"
[ -f $KIT/mpi_debug.patch ] || fail "kit/mpi_debug.patch missing (update the kit, see the handoff)"

# ---- 1. debug build (c1 source + prints) ----------------------------------------------------
if [ ! -x $M/mpi_c1_dbg ]; then
  cd $M/wt
  git checkout -q -- refine/prog/refine_run.f90
  git apply $KIT/mpi_debug.patch || fail "patch does not apply"
  cd $M/build
  echo "building debug binary ..."
  make -j$(nproc) > make_dbg.log 2>&1 || { git -C $M/wt checkout -q -- refine/prog/refine_run.f90; fail "make, see $M/build/make_dbg.log"; }
  cp suite/prog/discus_suite $M/mpi_c1_dbg
  git -C $M/wt checkout -q -- refine/prog/refine_run.f90   # leave the scratch source clean
fi

# ---- 2. one-cycle MPI fit with the debug binary -----------------------------------------------
T=$M/dbg_$(date +%H%M%S)
mkdir -p $T && cd $T
cp $KIT/refine/gauss.data .
cat > work_mpi.mac <<'EOF'
branch kuplot
func P_amp*exp(0.0-((r[0]-0.7)/1.3)**2) + P_bg, -5, 5, 0.1
exit
finished
EOF
cat > main.mac <<'EOF'
set prompt,off
seed 4711, 271828, 314159, 161803, 141421, 173205, 223606, 99991
refine
data xy, gauss.data
newparam P_amp, value:6.0, points:5
newparam P_bg,  value:0.2, points:5
set log, on
set cycle, 1
run work_mpi.mac
exit
exit
EOF
script -qec "timeout 400 mpirun ${MPIRUN_OPTS:-} --oversubscribe -n 9 $M/mpi_c1_dbg -macro main.mac < /dev/null" /dev/null > out.txt 2>&1
CLEAN() { sed 's/\x1b\[[0-9;]*m//g' "$@"; }

echo "## A. what the master receives (DBG lines)"
CLEAN out.txt | grep -E "DBG" | head -40
echo
echo "## B. errors"
CLEAN out.txt | grep -E "\*\*\*" | grep -v "Version\|\*\*\*\*\*\*" | head -10
echo
echo "## C. worker result files (sizes; they must differ from each other)"
ls -l DISCUS_SUITE_DERIVATIVES/data.* 2>&1 | awk '{print $5, $9}'
md5sum DISCUS_SUITE_DERIVATIVES/data.* 2>&1 | awk '{print $1}' | sort | uniq -c | sort -rn | head -3
echo
echo "## D. last lines of two worker logs"
for f in DISCUS_SUITE_DERIVATIVES/LOGFILE.0001 DISCUS_SUITE_DERIVATIVES/LOGFILE.0005; do
  echo "-- $f"; CLEAN $f | grep -v '^\s*$' | tail -n 6
done
echo
echo "(raw outputs: $T)"
echo "RESULT: DONE (information only)"

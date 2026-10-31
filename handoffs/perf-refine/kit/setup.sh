#!/bin/bash
# Builds discus_suite for the upstream base and for every commit of perf/refine.
# Run inside WSL:  bash ~/refine_validation/src/handoffs/perf-refine/kit/setup.sh
# Result: ~/refine_validation/bins/c0 .. c7 (+ cN.sha), kit copied to ~/refine_validation/kit
set -euo pipefail
ROOT=~/refine_validation
SRC=$ROOT/src
KITSRC=$SRC/handoffs/perf-refine/kit
cd $ROOT

mapfile -t SHAS < <(grep -v '^#' $KITSRC/COMMITS | awk '{print $1}')
[ ${#SHAS[@]} -eq 8 ] || { echo "COMMITS must list 8 commits (c0..c7)"; exit 1; }

rm -rf $ROOT/kit && cp -r $KITSRC $ROOT/kit
mkdir -p $ROOT/bins $ROOT/build

if [ ! -d $ROOT/wt ]; then
  git -C $SRC worktree add --detach $ROOT/wt ${SHAS[0]}
fi

export FC=$(which gfortran)
export PGPLOT_DIR=${PGPLOT_DIR:-/usr/lib/pgplot5}
cd $ROOT/build
if [ ! -f CMakeCache.txt ]; then
  cmake $ROOT/wt -DDIFFUSE_GRAPHICS=OFF -DCMAKE_BUILD_TYPE=Release > cmake.log 2>&1 \
    || { echo "cmake failed, see $ROOT/build/cmake.log"; exit 1; }
fi

for i in "${!SHAS[@]}"; do
  sha=${SHAS[$i]}
  if [ -x $ROOT/bins/c$i ] && [ "$(cat $ROOT/bins/c$i.sha 2>/dev/null)" = "$sha" ]; then
    echo "c$i $sha already built"; continue
  fi
  git -C $ROOT/wt checkout -q --detach $sha
  echo "building c$i $(git -C $ROOT/wt log -1 --format='%h %s')"
  make -j$(nproc) > make_c$i.log 2>&1 || { echo "build of c$i failed, see $ROOT/build/make_c$i.log"; exit 1; }
  cp suite/prog/discus_suite $ROOT/bins/c$i
  echo $sha > $ROOT/bins/c$i.sha
done
echo "SETUP DONE: $(ls $ROOT/bins | grep -c -v sha) binaries in $ROOT/bins"

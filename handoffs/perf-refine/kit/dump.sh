#!/bin/bash
# dump.sh N : condensed raw outputs of binaries c(N-1) and cN, for a FAIL report
N=$1; R=~/refine_validation/runs
for d in c$((N-1)) c$N; do
  for f in $R/$d/*/out.txt $R/$d/*/A.out $R/$d/*/B.out; do
    [ -f "$f" ] || continue
    echo "== $f"
    grep -E 'Cyc|^ +[0-9]+ |P_|Value of|[*][*][*]|No valid' "$f" | grep -v 'Version\|\*\*\*\*\*' | tail -n 25
  done
done 2>&1 | head -n 400

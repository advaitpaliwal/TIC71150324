#!/bin/sh
# usage: go.sh MODE  -> 10 runs, 4 at a time
H=$(cd "$(dirname "$0")" && pwd); PY=~/Scratch/planet-hunt/triceratops/.venv/bin/python
seq 0 9 | xargs -P 4 -I '{}' sh -c 'cd "$1"/out && OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "$2" -B -W ignore "$1"/run.py "$3" "$4" > "$1"/out/run_"$3"_"$4".log 2>&1' sh "$H" "$PY" "$1" '{}'

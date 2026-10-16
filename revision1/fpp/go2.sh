#!/bin/sh
# Remaining jobs (runner was killed after full_0-7): full 8,9 and sel 0-9, 6 at a time.
H=$(cd "$(dirname "$0")" && pwd); PY=~/Scratch/planet-hunt/triceratops/.venv/bin/python
printf 'full 8\nfull 9\nsel 0\nsel 1\nsel 2\nsel 3\nsel 4\nsel 5\nsel 6\nsel 7\nsel 8\nsel 9\n' | xargs -P 6 -L 1 sh -c 'cd "$0"/out && OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "$1" -B -W ignore "$0"/run.py "$2" "$3" > "$0"/out/run_"$2"_"$3".log 2>&1' "$H" "$PY"

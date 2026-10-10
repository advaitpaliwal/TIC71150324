#!/bin/sh
set -eu
# Ten joint runs; at most three simultaneous Python workers.
seq 0 9 | xargs -P 3 -I '{}' sh -c 'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 work/triceratops/.venv/bin/python -B work/deep/71150324/fpp/scripts/run_fpp.py 71150324 8.422187911756147 1443.4437448219462 "$1" > "work/deep/71150324/fpp/run_$1.log" 2>&1' sh '{}'
work/triceratops/.venv/bin/python -B -I work/deep/71150324/fpp/scripts/aggregate.py

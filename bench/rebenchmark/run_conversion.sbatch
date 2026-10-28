#!/bin/bash
# One shared decoder and queue for a pilot, coordinator, or multi-node workers.
set -euo pipefail
mode=$1
output=$2
checkout=$3
expected_commit=$4
python=$5
workers=${6:-1}
test "$(git -C "$checkout" rev-parse HEAD)" = "$expected_commit"
test -z "$(git -C "$checkout" status --porcelain --untracked-files=no)"
export PYTHONPATH="$checkout/src:$checkout/bench"
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

if [[ "$mode" == workers ]]; then
    srun --ntasks-per-node=1 bash "$0" worker "$output" "$checkout" "$expected_commit" "$python" "$workers"
    exit
fi

export TMPDIR
TMPDIR=$(mktemp -d /tmp/mappa-conversion.XXXXXX)
trap 'rm -rf "$TMPDIR"' EXIT
date -u
hostname
lscpu
"$python" -m pip freeze
"$python" -c 'import platform, sys; print(sys.version); print(platform.platform())'
if [[ "$mode" == pilot ]]; then
    cd "$checkout"
    "$python" -m pytest -q -p no:cacheprovider tests/test_saved_event_decoding.py tests/test_fragment_conversion.py tests/test_shared_execution.py
    "$python" -u -m rebenchmark.decoder_campaign coordinator "$output" &
    coordinator_pid=$!
    while [[ ! -f "$output/endpoint.json" ]]; do
        kill -0 "$coordinator_pid"
        sleep 1
    done
    "$python" -u -m rebenchmark.decoder_campaign worker "$output" --workers "$workers"
    wait "$coordinator_pid"
elif [[ "$mode" == coordinator ]]; then
    "$python" -u -m rebenchmark.decoder_campaign coordinator "$output"
elif [[ "$mode" == worker ]]; then
    "$python" -u -m rebenchmark.decoder_campaign worker "$output" --workers "$workers"
else
    echo "Unknown mode: $mode" >&2
    exit 2
fi
date -u

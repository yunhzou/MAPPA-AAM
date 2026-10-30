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
    srun --ntasks-per-node=1 bash "$checkout/bench/rebenchmark/run_conversion.sbatch" worker "$output" "$checkout" "$expected_commit" "$python" "$workers"
    exit
fi

export TMPDIR
# /tmp is tmpfs on NRT: a large queue and its backup would consume job RAM.
scratch_root=/var/tmp
case "$(findmnt -n -o FSTYPE -T "$scratch_root")" in
    tmpfs|ramfs|nfs|nfs4|lustre) echo "Scratch must be node-local disk: $scratch_root" >&2; exit 1 ;;
esac
TMPDIR=$(mktemp -d "$scratch_root/mappa-conversion.XXXXXX")
completed=false
cleanup() {
    if $completed; then rm -rf "$TMPDIR";
    else echo "Retained failed-job scratch: $TMPDIR" >&2; fi
}
trap cleanup EXIT
trap 'exit 143' TERM
trap 'exit 130' INT
mkdir -p "$output"
"$python" -c 'import json,os,socket,time; from pathlib import Path; p=Path(os.environ["TMPDIR"]); print(json.dumps(dict(host=socket.gethostname(),scratch=str(p),job=os.environ.get("SLURM_JOB_ID"),utc=time.time())))' > "$output/scratch-${SLURM_JOB_ID:-local}-$(hostname)-$$.json"
date -u
hostname
lscpu
"$python" -m pip freeze
"$python" -c 'import platform, sys; print(sys.version); print(platform.platform())'
case "$mode" in
    prepare) "$python" -u -m rebenchmark.inputs "$output" ;;
    prepare-queue|pilot-check|control|service|report|memory-stage)
        "$python" -u -m rebenchmark.launch "$mode" "$output" --workers "$workers" ;;
    worker) "$python" -u -m rebenchmark.execution worker "$output" --workers "$workers" ;;
    *) echo "Unknown mode: $mode" >&2; exit 2 ;;
esac

date -u
completed=true

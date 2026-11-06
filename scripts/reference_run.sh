#!/usr/bin/env bash
# Generate golden reference outputs with the original C++ hifiasm binary.
#
# Usage: scripts/reference_run.sh [threads]
#
# All artefacts land in reference/ (git-ignored). The script is intentionally
# single threaded by default so the output is reproducible.
set -euo pipefail

repository_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
reference_binary="${HIFIASM_REFERENCE_BINARY:-/home/keith/projects/hifiasm/hifiasm}"
input_reads="${HIFIASM_REFERENCE_INPUT:-/home/keith/Downloads/chr11-2M.fa.gz}"
thread_number="${1:-1}"
output_directory="${repository_root}/reference"
log_directory="${output_directory}/logs"

mkdir -p "${output_directory}" "${log_directory}"

run_profile() {
    local profile_name="$1"
    shift
    echo "[reference] profile=${profile_name} threads=${thread_number}"
    "${reference_binary}" \
        -t "${thread_number}" \
        -o "${output_directory}/${profile_name}" \
        "$@" \
        "${input_reads}" \
        >"${log_directory}/${profile_name}.stdout.log" \
        2>"${log_directory}/${profile_name}.stderr.log"
}

# Default HiFi-only assembly: produces *.bp.p_ctg.gfa, *.bp.p_utg.gfa, ...
run_profile default

# Same pipeline plus corrected reads and overlaps, used as intermediate gates.
run_profile write_paf_ec --write-paf --write-ec

sha256sum "${output_directory}"/*.gfa "${output_directory}"/*.bed \
    >"${output_directory}/checksums.sha256"

echo "[reference] checksums written to ${output_directory}/checksums.sha256"

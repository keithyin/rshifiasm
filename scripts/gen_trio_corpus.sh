#!/usr/bin/env bash
# Generate the synthetic trio corpus: split the chr11-2M reads into two
# pseudo-haplotypes by read-name checksum parity and write the two name
# lists consumed by `-3`/`-4`. Deterministic; data lands in reference/
# (git-ignored).
set -euo pipefail
repository_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
input_reads="${HIFIASM_REFERENCE_INPUT:-/home/keith/Downloads/chr11-2M.fa.gz}"
output_directory="${repository_root}/reference"
mkdir -p "${output_directory}"
zcat "${input_reads}" | awk '/^>/ {print substr($1, 2)}' > "${output_directory}/trio_all_names.txt"
: > "${output_directory}/trio_l1.txt"
: > "${output_directory}/trio_l2.txt"
while IFS= read -r name; do
    checksum="$(printf '%s' "${name}" | cksum)"
    if (( ${checksum%% *} % 2 == 1 )); then
        printf '%s\n' "${name}" >> "${output_directory}/trio_l1.txt"
    else
        printf '%s\n' "${name}" >> "${output_directory}/trio_l2.txt"
    fi
done < "${output_directory}/trio_all_names.txt"
echo "[trio-corpus] l1=$(wc -l < "${output_directory}/trio_l1.txt") l2=$(wc -l < "${output_directory}/trio_l2.txt")"

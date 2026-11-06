#!/usr/bin/env bash
# Compare the Rust implementation against the golden C++ reference outputs.
#
# Usage: scripts/compare_outputs.sh <reference_prefix> <candidate_prefix>
#
# Every artefact of the reference prefix must exist under the candidate prefix
# and be byte identical. On the first difference the script prints the byte
# offset together with a short context window to make the regression cheap to
# localise.
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "usage: $0 <reference_prefix> <candidate_prefix>" >&2
    exit 2
fi

reference_prefix="$1"
candidate_prefix="$2"
reference_directory="$(dirname "${reference_prefix}")"
reference_stem="$(basename "${reference_prefix}")"
candidate_stem="$(basename "${candidate_prefix}")"

differences=0
for reference_file in "${reference_directory}/${reference_stem}".*; do
    suffix="${reference_file#"${reference_directory}/${reference_stem}"}"
    candidate_file="$(dirname "${candidate_prefix}")/${candidate_stem}${suffix}"

    if [[ ! -f "${candidate_file}" ]]; then
        echo "[compare] MISSING ${candidate_file}"
        differences=$((differences + 1))
        continue
    fi

    if cmp --silent "${reference_file}" "${candidate_file}"; then
        echo "[compare] identical ${suffix}"
        continue
    fi

    differences=$((differences + 1))
    echo "[compare] DIFFERS   ${suffix}"
    cmp "${reference_file}" "${candidate_file}" || true
    echo "--- reference context ---"
    cmp -l "${reference_file}" "${candidate_file}" 2>/dev/null | head -n 5 || true
done

if [[ ${differences} -ne 0 ]]; then
    echo "[compare] ${differences} artefact(s) differ"
    exit 1
fi

echo "[compare] all artefacts are byte identical"

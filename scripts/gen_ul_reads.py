#!/usr/bin/env python3
"""Generate synthetic ultra-long reads for the `--ul` gate.

The reads are sampled from the contigs of a default assembly, optionally
reverse-complemented, and corrupted with ~5% substitutions and indels. This
exercises the UL alignment/rescue path end to end without needing a real
ultra-long dataset.

Usage: gen_ul_reads.py <assembly.gfa> <out.fq> [coverage]

The random seed is fixed so the gate input is reproducible.
"""

import random
import sys

COMPLEMENT = str.maketrans("ACGT", "TGCA")
SEED = 20261110
MIN_LENGTH = 40000
MAX_LENGTH = 80000


def main() -> None:
    if len(sys.argv) not in (3, 4):
        sys.exit(__doc__)
    gfa, output, coverage = sys.argv[1], sys.argv[2], float(sys.argv[3]) if len(sys.argv) == 4 else 20.0

    template = []
    with open(gfa) as handle:
        for line in handle:
            if line.startswith("S\t"):
                template.append(line.rstrip("\n").split("\t")[2])
    genome = "".join(template)
    if not genome:
        sys.exit(f"no segment lines in {gfa}")

    rng = random.Random(SEED)
    read_count = max(1, int(len(genome) * coverage / 60000))
    with open(output, "w") as handle:
        for index in range(read_count):
            length = min(rng.randint(MIN_LENGTH, MAX_LENGTH), len(genome))
            start = rng.randint(0, len(genome) - length)
            fragment = genome[start:start + length]
            if rng.random() < 0.5:
                fragment = fragment.translate(COMPLEMENT)[::-1]
            bases = []
            position = 0
            while position < len(fragment):
                roll = rng.random()
                if roll < 0.04:
                    bases.append(rng.choice("ACGT"))
                    position += 1
                elif roll < 0.05:
                    position += 1
                elif roll < 0.06:
                    bases.append(fragment[position])
                    bases.append(rng.choice("ACGT"))
                    position += 1
                else:
                    bases.append(fragment[position])
                    position += 1
            read = "".join(bases)
            handle.write(f"@ul_{index}\n{read}\n+\n{'I' * len(read)}\n")

    print(f"{read_count} reads from {len(genome)} bp", file=sys.stderr)


if __name__ == "__main__":
    main()

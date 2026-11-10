#!/usr/bin/env python3
"""Generate synthetic yak k-mer dumps for the trio-binning gate.

The reads named in the paternal/maternal lists are assigned to their
haplotype; every 31-mer of those reads is written to the matching dump with
the fixed count 5 (the default `mid_cnt`), so each read classifies as a
private haplotype read. Reads in neither list contribute no k-mer and end up
ambiguous, exactly like a real trio run where the dumps only carry the
parents' k-mers.

Usage: gen_yak_dump.py <in.fa[.gz]> <paternal.lst> <maternal.lst> <out_prefix>

Writes <out_prefix>.hap1.yak and <out_prefix>.hap2.yak in the format read by
`yak_ch_restore_core` (Trio.cpp): magic "YAK\\2", the u32 triple
(k, pre, counter_bits), then per bucket a u32 capacity, a u32 size and the
size u64 keys `(hash >> pre) << 10 | count`.
"""

import gzip
import struct
import sys

K = 31
PRE = 10
COUNTER_BITS = 10
MASK = (1 << (2 * K)) - 1
SHIFT = 2 * (K - 1)
WRITTEN_COUNT = 5

NT4 = [4] * 256
for base, code in ((b"A", 0), (b"C", 1), (b"G", 2), (b"T", 3),
                   (b"a", 0), (b"c", 1), (b"g", 2), (b"t", 3)):
    NT4[base[0]] = code


def yak_hash64(key, mask):
    key = (~key + (key << 21)) & mask
    key ^= key >> 24
    key = ((key + (key << 3)) + (key << 8)) & mask
    key ^= key >> 14
    key = ((key + (key << 2)) + (key << 4)) & mask
    key ^= key >> 28
    return (key + (key << 31)) & mask


def kmer_hashes(sequence):
    forward = 0
    reverse = 0
    run = 0
    hashes = set()
    for base in sequence:
        code = NT4[base]
        if code < 4:
            forward = ((forward << 2) | code) & MASK
            reverse = (reverse >> 2) | ((3 - code) << SHIFT)
            run += 1
            if run >= K:
                hashes.add(yak_hash64(min(forward, reverse), MASK))
        else:
            run = 0
            forward = reverse = 0
    return hashes


def read_fasta(path):
    opener = gzip.open if path.endswith(".gz") else open
    name = None
    chunks = []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).encode()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks).encode()


def read_names(path):
    names = set()
    with open(path, "rt") as handle:
        for line in handle:
            fields = line.split()
            if fields:
                names.add(fields[0])
    return names


def write_dump(path, hashes):
    buckets = [[] for _ in range(1 << PRE)]
    for value in hashes:
        buckets[value & ((1 << PRE) - 1)].append(
            ((value >> PRE) << COUNTER_BITS) | WRITTEN_COUNT)
    with open(path, "wb") as handle:
        handle.write(b"YAK\x02")
        handle.write(struct.pack("<3I", K, PRE, COUNTER_BITS))
        for keys in buckets:
            capacity = 4
            while capacity < max(4, len(keys) * 2):
                capacity <<= 1
            handle.write(struct.pack("<2I", capacity, len(keys)))
            for key in keys:
                handle.write(struct.pack("<Q", key))


def main():
    if len(sys.argv) != 5:
        print(__doc__, file=sys.stderr)
        return 2
    fasta, paternal_list, maternal_list, prefix = sys.argv[1:]
    paternal_names = read_names(paternal_list)
    maternal_names = read_names(maternal_list)

    paternal = set()
    maternal = set()
    for name, sequence in read_fasta(fasta):
        if name in paternal_names:
            paternal |= kmer_hashes(sequence)
        elif name in maternal_names:
            maternal |= kmer_hashes(sequence)

    write_dump(prefix + ".hap1.yak", paternal)
    write_dump(prefix + ".hap2.yak", maternal)
    print("paternal k-mers: {}, maternal k-mers: {}".format(len(paternal), len(maternal)),
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())

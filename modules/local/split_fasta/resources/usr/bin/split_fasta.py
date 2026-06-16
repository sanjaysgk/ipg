#!/usr/bin/env python3
"""Split a FASTA file into M chunks at protein boundaries.

Replicates get_fasta_offsets() and set_up_directories() from
FragPipe's msfragger_pep_split.py (lines 98-137).

Usage:
    split_fasta.py --fasta database.fasta --num_chunks 4 --outdir split_db
"""

import argparse
import mmap
import pathlib
import re
import shutil


def _array_split(lst, num_parts):
    """Split a list into num_parts roughly equal sublists (pure-Python replacement for numpy.array_split).

    Distributes remainder elements across the first chunks, matching numpy behaviour.

    Args:
        lst: The list to split.
        num_parts: Number of sublists to produce.

    Returns:
        List of sublists.
    """
    n = len(lst)
    q, r = divmod(n, num_parts)
    result = []
    idx = 0
    for i in range(num_parts):
        size = q + (1 if i < r else 0)
        result.append(lst[idx:idx + size])
        idx += size
    return result


def get_fasta_offsets(fasta_path: pathlib.Path, num_parts: int):
    """Get byte offsets for splitting FASTA at protein boundaries.

    Uses memory-mapped I/O to find all '\\n>' boundaries in the FASTA file,
    then splits the list of protein start positions into num_parts groups
    (same algorithm as FragPipe's numpy.array_split).

    Args:
        fasta_path: Path to the input FASTA file.
        num_parts: Number of chunks to split into.

    Returns:
        List of slice objects with (start, stop) byte offsets for each chunk.
    """
    with fasta_path.open("rb") as fo:
        mm = mmap.mmap(fo.fileno(), 0, access=mmap.ACCESS_READ)
        pos = [0] + [e.start() + 1 for e in re.compile(rb"\n>").finditer(mm)]
        startpos = [e[0] for e in _array_split(pos, num_parts)]
        return [slice(s, e) for s, e in zip(startpos, startpos[1:] + [len(mm)])]


def split_fasta(fasta_path: pathlib.Path, num_parts: int, outdir: pathlib.Path):
    """Split FASTA into num_parts chunks in outdir/{0,1,...}/original_name.fasta.

    Each chunk directory contains a copy of the FASTA with a subset of proteins,
    preserving the original filename. Chunks are split at protein boundaries
    so no protein sequence is truncated.

    Args:
        fasta_path: Path to the input FASTA file.
        num_parts: Number of chunks to create.
        outdir: Output directory (will be created; removed if exists).

    Returns:
        List of Path objects pointing to each chunk FASTA file.
    """
    offsets = get_fasta_offsets(fasta_path, num_parts)

    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True)

    chunk_paths = []
    with fasta_path.open("rb") as fo:
        for i, s in enumerate(offsets):
            chunk_dir = outdir / str(i)
            chunk_dir.mkdir()
            chunk_fasta = chunk_dir / fasta_path.name
            with chunk_fasta.open("wb") as f:
                mm = mmap.mmap(
                    fo.fileno(),
                    0 if s.stop is None else s.stop,
                    access=mmap.ACCESS_READ,
                )
                mm.seek(s.start)
                shutil.copyfileobj(mm, f)
            chunk_paths.append(chunk_fasta)

    return chunk_paths


def main():
    parser = argparse.ArgumentParser(
        description="Split FASTA file into chunks at protein boundaries"
    )
    parser.add_argument(
        "--fasta", required=True, type=pathlib.Path, help="Input FASTA file"
    )
    parser.add_argument(
        "--num_chunks", required=True, type=int, help="Number of chunks to split into"
    )
    parser.add_argument(
        "--outdir",
        default=pathlib.Path("split_db"),
        type=pathlib.Path,
        help="Output directory (default: split_db)",
    )
    args = parser.parse_args()

    if args.num_chunks < 1:
        parser.error("--num_chunks must be >= 1")

    if not args.fasta.exists():
        parser.error(f"FASTA file not found: {args.fasta}")

    chunks = split_fasta(args.fasta, args.num_chunks, args.outdir)
    for chunk in chunks:
        print(chunk)


if __name__ == "__main__":
    main()

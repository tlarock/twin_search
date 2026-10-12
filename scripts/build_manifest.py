#!/usr/bin/env python3
"""Checksum manifest for the reproducibility dataset.

One line per file: path, bytes, sha256. Mirrors the convention of the
published archive's MANIFEST.tsv, which is what a data deposit needs so a
downloader can prove they got what was published.

Per-directory MANIFEST.tsv files are skipped rather than listed: the
distance-matrices one carries an extra `values` column that the top-level
cannot, so both are kept and the top-level does not try to checksum a file
that changes whenever it is regenerated.
"""
import argparse
import hashlib
import os
import sys


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fin:
        for block in iter(lambda: fin.read(1 << 20), b''):
            h.update(block)
    return path, os.path.getsize(path), h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', required=True)
    ap.add_argument('--out', default=None, help='default: <root>/MANIFEST.tsv')
    ap.add_argument('--cores', type=int, default=4)
    args = ap.parse_args()
    root = os.path.abspath(args.root)
    out = args.out or os.path.join(root, 'MANIFEST.tsv')

    files = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = sorted(d for d in dirnames if not d.startswith('.'))
        for fn in sorted(filenames):
            if fn.startswith('.') or fn == 'MANIFEST.tsv':
                continue
            files.append(os.path.join(dirpath, fn))

    print(f"{len(files)} file(s) under {root}", flush=True)
    if args.cores > 1:
        from multiprocessing import Pool
        with Pool(args.cores) as p:
            rows = p.map(sha256, files, chunksize=16)
    else:
        rows = [sha256(f) for f in files]

    total = 0
    with open(out, 'w') as fh:
        fh.write("file\tbytes\tsha256\n")
        for path, size, digest in sorted(rows):
            total += size
            fh.write(f"{os.path.relpath(path, root)}\t{size}\t{digest}\n")
    print(f"{len(rows)} entries, {total / 2**30:.2f} GiB -> {out}")

    by_dir = {}
    for path, size, _ in rows:
        key = os.path.relpath(os.path.dirname(path), root)
        n, b = by_dir.get(key, (0, 0))
        by_dir[key] = (n + 1, b + size)
    print()
    for key in sorted(by_dir):
        n, b = by_dir[key]
        print(f"  {key:<22} {n:>5} files  {b / 2**30:>8.2f} GiB")
    return 0


if __name__ == '__main__':
    sys.exit(main())

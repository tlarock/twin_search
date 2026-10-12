#!/usr/bin/env python3
"""Acceptance test for the re-run distance matrices, and the archive manifest.

Two independent checks per file:

  1. TRANSFER. sha256 of the local .f32.zst against the sha256 the converter
     recorded on della, in the .meta. Proves the rsync was clean.
  2. SCIENCE. Decompressed values against the PUBLISHED .f32.zst.
     - hyperPD and jaccard must be IDENTICAL. Both were already shown to
       reproduce bit-for-bit under the newer xgi, so any difference here is a
       fault in this campaign, not a known change.
     - hyperNS is EXPECTED to differ: the published values use the 45-dim
       feature vector that upstream deleted, ours the 21-dim one from
       arXiv:2503.16959v2. The magnitude is reported, not asserted.

Compare DECOMPRESSED values, never compressed bytes: this campaign ran
`zstd -T6` and multithreaded zstd splits blocks differently from the
single-threaded published run, so identical content gives different archives.

Writes MANIFEST.tsv (name, bytes, sha256, values) for the data deposit.
"""
import argparse
import hashlib
import os
import subprocess
import sys

import numpy as np

CHUNK = 1 << 26                       # 64 MB of compressed output per read
DISTANCES = ('hyperNS', 'hyperPD', 'jaccard')
EXPECT_IDENTICAL = ('hyperPD', 'jaccard')


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fin:
        for block in iter(lambda: fin.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


class F32Stream:
    """Decompressed bytes of a .f32.zst, read in aligned blocks."""

    def __init__(self, path):
        self.proc = subprocess.Popen(['zstd', '-dcq', path], stdout=subprocess.PIPE)
        self.buf = b''

    def read(self, nbytes):
        out = [self.buf] if self.buf else []
        have = len(self.buf)
        self.buf = b''
        while have < nbytes:
            block = self.proc.stdout.read(nbytes - have)
            if not block:
                break
            out.append(block)
            have += len(block)
        data = b''.join(out)
        # never split a float across a block boundary
        keep = len(data) - (len(data) % 4)
        self.buf = data[keep:]
        return data[:keep]

    def close(self):
        try:
            self.proc.stdout.close()
        finally:
            return self.proc.wait()


def signature(path):
    """Order-INSENSITIVE fingerprint of a matrix's values.

    Positional comparison against the published archive is meaningless here,
    and this is the fifth time that trap has been sprung in this campaign.
    BOTH orderings are thread artefacts:

      - within a projection, the twin order the search emitted, so the pair
        order in the condensed vector differs;
      - across projections, the ROW order of the enumeration itself. Measured
        2026-10-12: our m=5 and m=8 total_twins sequences differ from the
        published ones element-wise but are identical as multisets.

    What IS invariant is the multiset of values: the twin sets were verified
    identical as sets, so the set of unordered pairs is the same, so the bag
    of distances must be too.

    Returns (count, min, max, float64 sum, multiset hash). The hash sums a
    64-bit mix of each value's raw bit pattern mod 2^64 - commutative, so it
    is independent of order, and exact for multiset equality up to a collision
    probability far below anything else that could go wrong here. It needs one
    streaming pass and constant memory, which sorting 8.6 billion floats does
    not.
    """
    st = F32Stream(path)
    n = 0
    mn, mx, total, h = float('inf'), float('-inf'), 0.0, 0
    try:
        while True:
            block = st.read(CHUNK)
            if not block:
                break
            u = np.frombuffer(block, dtype='<u4')
            f = u.view('<f4')
            n += f.size
            mn = min(mn, float(f.min()))
            mx = max(mx, float(f.max()))
            total += float(f.astype(np.float64).sum())
            with np.errstate(over='ignore'):
                x = u.astype(np.uint64)
                x ^= x << np.uint64(33)
                x *= np.uint64(0xff51afd7ed558ccd)
                x ^= x >> np.uint64(29)
                x *= np.uint64(0xc4ceb9fe1a85ec53)
                x ^= x >> np.uint64(32)
                h = (h + int(x.sum(dtype=np.uint64))) & 0xFFFFFFFFFFFFFFFF
    finally:
        rc = st.close()
    if rc != 0:
        raise RuntimeError(f"zstd -d exited {rc} on {path}")
    return n, mn, mx, total, h


def check_one(name, args):
    ours = os.path.join(args.ours, name)
    pub = os.path.join(args.published, name)
    dist = name.rsplit('_', 1)[1].split('.')[0]

    size = os.path.getsize(ours)
    digest = sha256(ours)
    meta = {}
    if os.path.exists(ours + '.meta'):
        with open(ours + '.meta') as fh:
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                if len(parts) == 2:
                    meta[parts[0]] = parts[1]
    sha_ok = (meta.get('sha256') == digest) if 'sha256' in meta else None

    if not os.path.exists(pub):
        return dict(name=name, dist=dist, size=size, sha=digest, sha_ok=sha_ok,
                    values=int(meta.get('values', -1)), verdict='no published counterpart')

    na, amn, amx, asum, ah = signature(ours)
    nb, bmn, bmx, bsum, bh = signature(pub)
    if na != nb:
        verdict = f'LENGTH {na:,} vs published {nb:,}'
    elif ah == bh:
        verdict = 'identical multiset'
    else:
        # The sums are float64 over differently ordered float32 values, so
        # they disagree in the last ulps even for equal multisets; quote the
        # means, which is what the figures plot.
        verdict = (f'DIFFERS  mean {asum / na:.6f} vs {bsum / nb:.6f}  '
                   f'range [{amn:.4f},{amx:.4f}] vs [{bmn:.4f},{bmx:.4f}]')
    return dict(name=name, dist=dist, size=size, sha=digest, sha_ok=sha_ok,
                values=na, verdict=verdict)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ours', required=True)
    ap.add_argument('--published', required=True)
    ap.add_argument('--manifest', default=None)
    ap.add_argument('--distances', nargs='*', default=list(DISTANCES))
    ap.add_argument('--m-max', type=int, default=13)
    ap.add_argument('--cores', type=int, default=3)
    args = ap.parse_args()

    names = sorted(f for f in os.listdir(args.ours)
                   if f.endswith('.f32.zst')
                   and f.rsplit('_', 1)[1].split('.')[0] in args.distances
                   and int(f.split('_m-')[1].split('_')[0]) <= args.m_max)
    print(f"{len(names)} file(s)\n")

    if args.cores > 1:
        from multiprocessing import Pool
        with Pool(args.cores) as p:
            rows = p.starmap(check_one, [(n, args) for n in names])
    else:
        rows = [check_one(n, args) for n in names]

    rows.sort(key=lambda r: (r['dist'], int(r['name'].split('_m-')[1].split('_')[0])))
    print(f"  {'m':<4} {'distance':<9} {'values':>14}  {'sha':<5}  verdict")
    bad = 0
    for r in rows:
        m = r['name'].split('_m-')[1].split('_')[0]
        shatxt = 'ok' if r['sha_ok'] else ('--' if r['sha_ok'] is None else 'BAD')
        flag = ''
        if r['dist'] in EXPECT_IDENTICAL and r['verdict'] != 'identical multiset':
            flag = '   <-- EXPECTED IDENTICAL'
            bad += 1
        if r['sha_ok'] is False:
            flag += '   <-- SHA MISMATCH'
            bad += 1
        print(f"  {m:<4} {r['dist']:<9} {r['values']:>14,}  {shatxt:<5}  {r['verdict']}{flag}")

    if args.manifest:
        with open(args.manifest, 'w') as fh:
            fh.write("file\tbytes\tsha256\tvalues\n")
            for r in sorted(rows, key=lambda x: x['name']):
                fh.write(f"{r['name']}\t{r['size']}\t{r['sha']}\t{r['values']}\n")
        print(f"\nmanifest -> {args.manifest}")

    print(f"\n{'FAILURES: ' + str(bad) if bad else 'all checks passed'}")
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())

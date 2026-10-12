#!/usr/bin/env python3
"""Convert ONE gzipped decimal distance matrix to flat float32 + zstd.

Adapted from gram_mates results/python/validation/convert_distance_matrices.py,
which reads members out of a 205 GB tar. Ours are loose files, one per
(m, distance), so this takes a single path and is driven by a SLURM array -
the tar version had to be serial because gzip cannot be seeked.

Format (unchanged, so the published readers work): a flat little-endian
float32 array, NO row delimiters and NO index. Row r is the next
total_twins[r] * (total_twins[r] - 1) / 2 values, which every reader already
knows from the exhaustive CSV. float32 is safe because the figures plot
per-projection mean and std, which move by <= 8.5e-09; uint8 measured ~4e-4
and was rejected as visible at plot precision.

zstd level is per distance - 19 for hyperNS/hyperPD, 12 for jaccard - which is
an ENCODER setting only: same frame, same reader, same `zstd -d`. jaccard is
234M values from ~150 distinct ones, and long-range matching at 19 costs 15x
the time for ~21% (measured at m=10: 179 s vs 12 s).
"""
import argparse
import hashlib
import os
import re
import resource
import shutil
import subprocess
import sys
import time
import warnings

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'analysis', 'python'))
from distance_stats import open_text  # noqa: E402

CHUNK = 1 << 24
LEVELS = {'hyperNS': 19, 'hyperPD': 19, 'jaccard': 12}
DEFAULT_LEVEL = 19
SUFFIX = '.csv.gz'
EDGE = 10_000          # values compared at each end during verification
NAME_RE = re.compile(r'^n-(\d+)_m-(\d+)_k-(\d+)_(.*)_(hyperNS|hyperPD|jaccard)$')


# ---------------------------------------------------------------- parsing ---
def _make_parser():
    """Bind the text->float parser once, preferring the fast path.

    np.fromstring(sep=',') is the only vectorised text parser numpy exposes,
    but its BINARY mode is deprecated and the whole function has been on the
    way out for years. Probe it with warnings-as-errors so a numpy that has
    finally removed it degrades to a slower parser instead of killing a job
    that is an hour in. della runs numpy 2.4.6; this is checked there, not
    assumed from the laptop's 2.2.6.
    """
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            probe = np.fromstring("1.5,2.5,3.5", sep=",")
        if probe.size == 3 and probe[0] == 1.5 and probe[2] == 3.5:
            return ("np.fromstring", lambda s: np.fromstring(s, sep=","))
    except Exception:
        pass
    return ("np.array(split)", lambda s: np.array(s.split(","), dtype=np.float64))


PARSER_NAME, parse_values = _make_parser()


# ------------------------------------------------------------------ zstd ----
def pick_zstd(staged=None):
    """Resolve a WORKING zstd, same tiers as slurm-scripts/zstd_tool.sh.

    della has /usr/bin/zstd on the LOGIN node but not on compute nodes, which
    is what killed the first exact_mates array (18 tasks, 00:00:00, "zstd:
    command not found"). The staged copy on the shared filesystem is the tier
    that actually fires in a batch job.
    """
    if shutil.which("zstd"):
        return ["zstd"]
    staged = staged or os.path.join(HERE, '..', 'slurm-scripts', 'bin', 'zstd')
    staged = os.path.abspath(staged)
    if os.access(staged, os.X_OK):
        try:
            subprocess.run([staged, "-q", "-c"], input=b"x",
                           stdout=subprocess.DEVNULL, check=True)
            return [staged]
        except Exception:
            pass
    sys.exit(f"no working zstd on {os.uname().nodename}: tried PATH and {staged}")


# ------------------------------------------------------------ exhaustive ----
def total_twins_only(path):
    """total_twins per projection, without parsing the twin hypergraphs.

    Each row is header / sizedist / "unfilt,filt" / ... / twins, so summing the
    unfiltered counts over the middle fields is the whole job. read_data would
    also build 1.2M hypergraphs at m=13 to recover one integer per row.
    """
    out = []
    with open_text(path) as fin:
        for line in fin:
            parts = line.rstrip("\n").split("/")
            total = 0
            for i in range(2, len(parts) - 1, 2):
                total += int(parts[i].split(",")[0])
            out.append(total)
    return out


def classify(got, lengths):
    """-> (kind, detail) for a value count that missed the prediction.

    A file can legitimately hold FEWER projections than the exhaustive CSV if
    the computation was cut short: the count then matches a PREFIX of the row
    lengths. That is how several PUBLISHED files look. Ours are complete runs,
    so 'short' here means a task died rather than a known property of the
    data - the caller refuses to publish either way.
    """
    run = 0
    for i, n in enumerate(lengths):
        if run == got:
            return ("short", f"holds {i} of {len(lengths)} projections; "
                             f"the last {len(lengths) - i} were never written")
        run += n
    if run == got:
        return ("ok", "")
    return ("misaligned", f"count {got:,} matches no prefix of the row lengths")


# ------------------------------------------------------------- conversion ---
def convert(src, dst, zcmd, jobs, level):
    """Stream gzipped decimal text -> flat float32 -> zstd.

    Returns (values_written, first_edge, last_edge). The edges are kept so the
    verify pass can prove the OUTPUT holds them: the count alone is counted on
    the INPUT side of the pipe and would not notice zstd truncating.
    """
    import gzip
    out = open(dst, 'wb')
    proc = subprocess.Popen(zcmd + [f'-{level}', f'-T{jobs}', '-q', '-c'],
                            stdin=subprocess.PIPE, stdout=out)
    total = 0
    tail = b''
    first = np.empty(0, dtype='<f4')
    last = np.empty(0, dtype='<f4')

    def emit(a, body_len):
        nonlocal total, first, last
        # np.fromstring stops silently at the first value it cannot read, so a
        # malformed byte would truncate the row and misalign every row after
        # it. Count the separators and insist the parse saw them all.
        if a.size != body_len:
            raise ValueError(f"parsed {a.size} values from a block holding "
                             f"{body_len}: the source is malformed")
        f4 = a.astype('<f4')
        proc.stdin.write(f4.tobytes())
        if first.size < EDGE:
            first = np.concatenate([first, f4[:EDGE - first.size]])
        last = f4[-EDGE:] if f4.size >= EDGE else np.concatenate([last, f4])[-EDGE:]
        total += a.size

    try:
        with gzip.open(src, 'rb') as fin:
            while True:
                chunk = fin.read(CHUNK)
                if not chunk:
                    break
                buf = tail + chunk
                # Split only at a separator, so a number straddling the chunk
                # boundary is carried whole into the next round.
                cut = max(buf.rfind(b','), buf.rfind(b'\n'))
                if cut < 0:
                    tail = buf
                    continue
                tail = buf[cut + 1:]
                body = buf[:cut].replace(b'\n', b',').strip(b',')
                if not body:
                    continue
                text = body.decode()
                emit(parse_values(text), text.count(',') + 1)
        rest = tail.strip().strip(b',')
        if rest:
            text = rest.decode()
            emit(parse_values(text), text.count(',') + 1)
    finally:
        proc.stdin.close()
        rc = proc.wait()
        out.close()
    if rc != 0:
        raise RuntimeError(f"zstd exited {rc}; {dst} is truncated")
    return total, first, last


def verify(path, zcmd, want_values, first, last):
    """Read the COMPRESSED file back: byte count, both edges, sha256.

    The conversion counts what it fed INTO the pipe. Decompression is where a
    truncated or corrupt frame shows up, and it is nearly free (measured
    across this campaign), so there is no reason to skip it.
    """
    h = hashlib.sha256()
    with open(path, 'rb') as fin:
        for block in iter(lambda: fin.read(1 << 20), b''):
            h.update(block)
    proc = subprocess.Popen(zcmd + ['-dcq', path], stdout=subprocess.PIPE)
    nbytes = 0
    head = b''
    tailbuf = b''
    try:
        for block in iter(lambda: proc.stdout.read(1 << 20), b''):
            nbytes += len(block)
            if len(head) < EDGE * 4:
                head += block[:EDGE * 4 - len(head)]
            tailbuf = (tailbuf + block)[-EDGE * 4:]
    finally:
        proc.stdout.close()
        rc = proc.wait()
    if rc != 0:
        raise RuntimeError(f"zstd -d exited {rc}: {path} does not decompress")
    if nbytes != want_values * 4:
        raise RuntimeError(f"{path} holds {nbytes // 4:,} values, expected "
                           f"{want_values:,}")
    got_first = np.frombuffer(head, dtype='<f4')
    got_last = np.frombuffer(tailbuf, dtype='<f4')
    if not np.array_equal(got_first, first[:got_first.size]):
        raise RuntimeError(f"{path}: leading values differ from the source")
    if not np.array_equal(got_last, last[-got_last.size:]):
        raise RuntimeError(f"{path}: trailing values differ from the source")
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--src', required=True, help='one <base>_<dist>.csv.gz')
    ap.add_argument('--dest', required=True, help='output directory')
    ap.add_argument('--exhaustive-dir', required=True)
    ap.add_argument('--jobs', type=int, default=6, help='zstd threads')
    ap.add_argument('--level', type=int, default=None, help='override LEVELS')
    ap.add_argument('--force', action='store_true')
    args = ap.parse_args()

    src = os.path.abspath(args.src)
    base = os.path.basename(src)
    if not base.endswith(SUFFIX):
        sys.exit(f"expected a {SUFFIX} file, got {base}")
    stem = base[:-len(SUFFIX)]
    mo = NAME_RE.match(stem)
    if not mo:
        sys.exit(f"cannot read n/m/k/distance out of {stem}")
    n, m, k, tag, dist = int(mo.group(1)), int(mo.group(2)), int(mo.group(3)), \
        mo.group(4), mo.group(5)
    level = args.level if args.level is not None else LEVELS.get(dist, DEFAULT_LEVEL)

    os.makedirs(args.dest, exist_ok=True)
    dst = os.path.join(args.dest, stem + '.f32.zst')
    if os.path.exists(dst) and not args.force:
        print(f"exists, skipping: {dst}  (use --force)")
        return 0

    zcmd = pick_zstd()
    exh = os.path.join(args.exhaustive_dir, f"n-{n}_m-{m}_k-{k}_{tag}.csv")
    print(f"=== {stem}  level -{level}  zstd -T{args.jobs}  parser {PARSER_NAME} ===",
          flush=True)
    print(f"    zstd: {zcmd[0]}", flush=True)

    tw = total_twins_only(exh)
    lengths = [t * (t - 1) // 2 for t in tw]
    expect = sum(lengths)
    print(f"    {len(tw)} projections, expecting {expect:,} values", flush=True)

    # Write to a temp name and rename only once the file has been read back.
    # A half-written .f32.zst is indistinguishable from a good one by name,
    # and every downstream reader trusts the name.
    tmp = dst + f'.tmp.{os.getpid()}'
    t0 = time.time()
    try:
        total, first, last = convert(src, tmp, zcmd, args.jobs, level)
        kind, detail = ("ok", "") if total == expect else classify(total, lengths)
        if kind != "ok":
            raise RuntimeError(f"{kind}: {detail}")
        sha = verify(tmp, zcmd, total, first, last)
        os.replace(tmp, dst)
    except BaseException as exc:
        if os.path.exists(tmp):
            os.remove(tmp)
        print(f"=== FAILED {stem}: {exc} ===", file=sys.stderr, flush=True)
        return 1
    secs = int(time.time() - t0)

    src_sz, dst_sz = os.path.getsize(src), os.path.getsize(dst)
    rss = max(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
              resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
    if sys.platform == 'darwin':      # macOS reports bytes, Linux kilobytes
        rss //= 1024

    with open(dst + '.meta', 'w') as fh:
        for key, val in [
            ('status', 'ok'), ('kind', 'f32-convert'), ('distance', dist),
            ('n', n), ('m', m), ('k', k),
            ('projections', len(tw)), ('values', total),
            ('src_bytes', src_sz), ('dst_bytes', dst_sz),
            ('ratio_vs_gz', f"{src_sz / max(dst_sz, 1):.2f}"),
            ('level', level), ('zstd_threads', args.jobs), ('zstd_cmd', zcmd[0]),
            ('parser', PARSER_NAME), ('numpy', np.__version__),
            ('seconds', secs), ('peak_rss_kb', rss),
            ('sha256', sha),
            ('node', os.uname().nodename),
            ('job', os.environ.get('SLURM_JOB_ID', '-')),
            ('source', base),
        ]:
            fh.write(f"{key}\t{val}\n")

    print(f"    {src_sz / 1e9:.2f} GB gz -> {dst_sz / 1e9:.2f} GB f32.zst "
          f"({src_sz / max(dst_sz, 1):.2f}x)  {total:,} values  "
          f"{secs}s  peak {rss / 1048576:.1f} GB", flush=True)
    print(f"    sha256 {sha}", flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())

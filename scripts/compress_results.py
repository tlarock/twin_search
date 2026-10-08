#!/usr/bin/env python3
"""Compress finished result cells with zstd, in place and verified.

The sampled dataset is ~29 GB of CSV that compresses ~38x, because a twin list
is a few thousand near-identical integer tuples. That ratio is the difference
between a dataset that fits in a backup and one that does not: /scratch on the
cluster is not backed up, and the plain form is too large to move over a
domestic link.

Readers do not need to change. `distance_stats.resolve_path` answers a request
for "<cell>.csv" with "<cell>.csv.zst" when that is what exists, so notebooks
and analysis scripts keep naming the plain file.

WHAT THIS WILL NOT COMPRESS, and why it matters more than what it will:

  - *.partial-*rows - a censored cell, and the input to its own resume. The
    resume machinery appends to plain text and must never meet a .zst.
  - anything under a .resume-* staging directory, for the same reason.
  - *.lock, *.superseded, *.tmp - transient or already-dead.
  - a cell with fewer rows than its own "samples-N" filename claims. A short
    cell is one that is still being worked on, or one that failed; either way
    compressing it would freeze a partial result into something that looks
    finished.
  - a cell whose .meta says anything other than "status ok".

The original is removed only after the archive has been decompressed again and
compared byte for byte with it. zstd's own --rm does not do this: it trusts the
write. For a dataset that is about to become the only copy, a verify pass at
~2 GB/s is worth the few minutes it costs.
"""

import argparse
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

# "n-9_m-43_k-3_samples-500.csv" -> 500 rows expected.
SAMPLES_RE = re.compile(r"_samples-(\d+)")

SKIP_SUBSTRINGS = (".partial-", ".lock", ".superseded", ".tmp")


def classify(path):
    """-> (ok_to_compress, reason). Reason is printed when ok is False."""
    base = os.path.basename(path)

    for bad in SKIP_SUBSTRINGS:
        if bad in base:
            return False, f"transient or censored ({bad})"
    if os.sep + ".resume-" in path:
        return False, "inside a resume staging directory"

    meta = path + ".meta"
    if os.path.exists(meta):
        with open(meta) as fh:
            status = dict(
                line.rstrip("\n").split("\t", 1)
                for line in fh if "\t" in line
            ).get("status", "?")
        if status != "ok":
            return False, f"meta says status={status}"

    m = SAMPLES_RE.search(base)
    if m:
        want = int(m.group(1))
        have = count_rows(path)
        if have < want:
            return False, f"short: {have} rows < {want} the name claims"
        return True, ""

    # Exhaustive enumerations carry no row target - completeness is the exit
    # code, which only the .meta records. With no meta there is nothing here
    # that can distinguish a finished file from an interrupted one.
    if os.path.exists(meta):
        return True, ""
    return False, "no samples-N in the name and no meta: completeness unknown"


def show(path):
    """Shortest readable form: relative to cwd when that is not a walk upwards."""
    rel = os.path.relpath(path)
    return path if rel.startswith("..") else rel


def count_rows(path):
    with open(path, "rb") as fh:
        return sum(buf.count(b"\n") for buf in iter(lambda: fh.read(1 << 20), b""))


def compress(path, level, keep):
    """Compress one file, verifying before the original is removed."""
    tmp = path + ".zst.tmp"
    final = path + ".zst"
    before = os.path.getsize(path)
    try:
        subprocess.run(["zstd", f"-{level}", "-q", "-T1", "-o", tmp, "-f", path],
                       check=True)
        # Decompress the archive and diff it against the file it came from.
        # cmp reads both streams; nothing is held in memory.
        with subprocess.Popen(["zstd", "-dcq", tmp],
                              stdout=subprocess.PIPE) as proc:
            rc = subprocess.run(["cmp", "-s", "-", path],
                                stdin=proc.stdout).returncode
            proc.stdout.close()
            if proc.wait() not in (0, -13):
                raise RuntimeError("decompression failed during verify")
        if rc != 0:
            raise RuntimeError("round trip differs from the original")
        os.replace(tmp, final)
        after = os.path.getsize(final)
        if not keep:
            os.remove(path)
        return before, after, None
    except Exception as exc:                      # noqa: BLE001 - reported, not raised
        if os.path.exists(tmp):
            os.remove(tmp)
        return before, 0, f"{type(exc).__name__}: {exc}"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dirs", nargs="+", help="results directories to walk")
    ap.add_argument("-l", "--level", type=int, default=19,
                    help="zstd level (default 19, matching the published "
                         "dataset; -22 --ultra measured WORSE on this data)")
    ap.add_argument("-j", "--jobs", type=int, default=4,
                    help="files compressed at once (default 4)")
    ap.add_argument("-n", "--dry-run", action="store_true")
    ap.add_argument("--keep", action="store_true",
                    help="leave the plain file in place as well")
    args = ap.parse_args()

    todo, skipped = [], []
    for root_dir in args.dirs:
        for root, _, files in os.walk(root_dir):
            for name in sorted(files):
                if not name.endswith(".csv"):
                    continue
                path = os.path.join(root, name)
                if os.path.exists(path + ".zst"):
                    skipped.append((path, "already has a .zst beside it"))
                    continue
                ok, why = classify(path)
                (todo if ok else skipped).append(path if ok else (path, why))

    if skipped:
        print(f"skipping {len(skipped)} file(s):")
        for path, why in sorted(skipped, key=lambda s: s[1]):
            print(f"  {why:<48} {show(path)}")
        print()

    if not todo:
        print("nothing to compress")
        return 0

    # Largest first: with a 1.2 GB cell in the set and ~6 minutes of
    # single-threaded work in it, starting it last would leave every other
    # worker idle waiting for it.
    todo.sort(key=os.path.getsize, reverse=True)
    total = sum(os.path.getsize(p) for p in todo)
    print(f"{len(todo)} file(s), {total / 2**30:.2f} GiB, "
          f"zstd -{args.level} on {args.jobs} job(s)")
    if args.dry_run:
        for path in todo[:10]:
            print(f"  {os.path.getsize(path) / 2**20:9.1f} MiB  "
                  f"{show(path)}")
        if len(todo) > 10:
            print(f"  ... and {len(todo) - 10} more")
        return 0

    done = before_sum = after_sum = 0
    failures = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for path, (before, after, err) in zip(todo, pool.map(
                lambda p: compress(p, args.level, args.keep), todo)):
            done += 1
            if err:
                failures.append((path, err))
                print(f"  [{done}/{len(todo)}] FAILED {show(path)}: {err}",
                      file=sys.stderr)
                continue
            before_sum += before
            after_sum += after
            print(f"  [{done}/{len(todo)}] {before / 2**20:8.1f} -> "
                  f"{after / 2**20:7.2f} MiB  {before / max(after, 1):6.1f}x  "
                  f"{show(path)}", flush=True)

    print(f"\n{done - len(failures)} compressed: "
          f"{before_sum / 2**30:.2f} -> {after_sum / 2**30:.3f} GiB "
          f"({before_sum / max(after_sum, 1):.1f}x overall)")
    if failures:
        print(f"{len(failures)} FAILED - originals left in place:", file=sys.stderr)
        for path, err in failures:
            print(f"  {show(path)}: {err}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

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


EXHAUSTIVE_RE = re.compile(
    r"^n-(\d+)_m-(\d+)_k-(\d+)_non-uniform_exhaustive_projections\.csv$")

_EXPECTED = None


def exhaustive_expected(basename):
    """-> rows an exhaustive enumeration must have, or None if not tabulated.

    An exhaustive run's row count is a DETERMINISTIC function of (n, m, k) -
    the number of non-isomorphic projections. That makes it exactly checkable,
    unlike a sampled cell where only a lower bound is meaningful.
    """
    global _EXPECTED
    if _EXPECTED is None:
        _EXPECTED = {}
        table = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "exhaustive_row_counts.tsv")
        try:
            with open(table) as fh:
                for line in fh:
                    if line.startswith("#") or line.startswith("n\t"):
                        continue
                    n, m, k, rows = line.split()
                    _EXPECTED[(int(n), int(m), int(k))] = int(rows)
        except OSError:
            pass
    mo = EXHAUSTIVE_RE.match(basename)
    if not mo:
        return None
    n, m, k = (int(g) for g in mo.groups())
    return _EXPECTED.get((n, m, k))


def expected_rows(basename):
    """-> (expected, exact) from the FILENAME alone, or (None, _).

    `exact` says whether the count must match rather than merely be reached:
    a sampled cell can legitimately hold more rows than its target, an
    exhaustive enumeration cannot.
    """
    mo = SAMPLES_RE.search(basename)
    if mo:
        return int(mo.group(1)), False
    want = exhaustive_expected(basename)
    return (want, True) if want is not None else (None, True)


def classify(path):
    """-> (ok_to_compress, reason). Reason is printed when ok is False.

    THE DATA OUTRANKS THE METADATA. A row count read off the file is evidence;
    a .meta is a claim written by a signal handler under time pressure, and it
    can be - has been - wrong about a file that is perfectly good.

    The case that taught this: an exhaustive cell whose SIGTERM trap stamped
    `status timeout-empty, rows 0` at 07:55 while the child was still running.
    The child finished at 13:15, so the file's mtime was five hours LATER than
    the metadata declaring it empty. The cell was complete and byte-identical
    to the published dataset, and an earlier version of this function refused
    to compress it on the meta's say-so - and was very nearly taken as grounds
    to DELETE 404 MB of correct results.
    """
    base = os.path.basename(path)

    # Name-based refusals come first. These are about what the file IS, not
    # what it contains: a partial is an input to its own resume whether or not
    # it happens to hold a full complement of rows.
    for bad in SKIP_SUBSTRINGS:
        if bad in base:
            return False, f"transient or censored ({bad})"
    if os.sep + ".resume-" in path:
        return False, "inside a resume staging directory"

    meta_status = None
    if os.path.exists(path + ".meta"):
        with open(path + ".meta") as fh:
            meta_status = dict(
                line.rstrip("\n").split("\t", 1)
                for line in fh if "\t" in line
            ).get("status", "?")

    want, exact = expected_rows(base)
    if want is not None:
        have = count_rows(path)
        if exact and have != want:
            return False, f"{have} rows, expected exactly {want}"
        if not exact and have < want:
            return False, f"short: {have} rows < {want} the name claims"
        if meta_status not in (None, "ok"):
            # Deliberately NOT a refusal. The file passed the only check that
            # looks at content; say the meta disagrees and move on.
            print(f"  note: {show(path)} passes its row check ({have}) but its "
                  f"meta says status={meta_status} - STALE META, compressing anyway")
        return True, ""

    # No check derivable from the name. Only now does the meta get a vote.
    if meta_status == "ok":
        return True, ""
    if meta_status is not None:
        return False, f"no row check available and meta says status={meta_status}"
    return False, "no expected row count and no meta: completeness unknown"


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


def decompress(path):
    """The inverse of compress(), same verify-before-delete discipline."""
    plain = path[:-len(".zst")]
    tmp = plain + ".tmp"
    try:
        subprocess.run(["zstd", "-dq", "-o", tmp, "-f", path], check=True)
        os.replace(tmp, plain)
        os.remove(path)
        return os.path.getsize(plain), None
    except Exception as exc:                      # noqa: BLE001
        if os.path.exists(tmp):
            os.remove(tmp)
        return 0, f"{type(exc).__name__}: {exc}"


def verify(path):
    """Decompress a cell and check it still holds the rows its name claims.

    This is the check worth running before the dataset becomes the only copy.
    A .zst that decodes cleanly can still be the WRONG cell - a half-finished
    compression renamed by hand, a file copied over another - and the row
    count catches that where a checksum of the archive against itself cannot.
    Exhaustive cells are checked for EXACT equality against
    exhaustive_row_counts.tsv; sampled cells only for a lower bound.
    """
    want, exact = expected_rows(os.path.basename(path)[:-len(".zst")])
    try:
        with subprocess.Popen(["zstd", "-dcq", path],
                              stdout=subprocess.PIPE) as proc:
            rows = sum(buf.count(b"\n")
                       for buf in iter(lambda: proc.stdout.read(1 << 20), b""))
            if proc.wait() != 0:
                return rows, "archive did not decompress cleanly"
    except Exception as exc:                      # noqa: BLE001
        return 0, f"{type(exc).__name__}: {exc}"
    if want is None:
        return rows, None
    if exact and rows != want:
        return rows, f"{rows} rows, expected exactly {want}"
    if not exact and rows < want:
        return rows, f"{rows} rows < {want} the name claims"
    return rows, None


def run_verify(dirs, jobs):
    todo = [os.path.join(root, name)
            for d in dirs for root, _, files in os.walk(d)
            for name in sorted(files) if name.endswith(".csv.zst")]
    if not todo:
        print("no compressed cells found")
        return 0
    print(f"verifying {len(todo)} compressed cell(s) on {jobs} job(s)")
    bad, total_rows = [], 0
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        for path, (rows, err) in zip(todo, pool.map(verify, todo)):
            total_rows += rows
            if err:
                bad.append((path, err))
                print(f"  BAD  {show(path)}: {err}", file=sys.stderr)
    print(f"{len(todo) - len(bad)}/{len(todo)} cells intact, "
          f"{total_rows:,} rows total")
    if bad:
        print(f"{len(bad)} FAILED verification", file=sys.stderr)
        return 1
    return 0


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
    ap.add_argument("--verify", action="store_true",
                    help="decompress every .csv.zst and check its row count "
                         "against the samples-N in its name; change nothing")
    ap.add_argument("--decompress", action="store_true",
                    help="the inverse: .csv.zst back to .csv")
    args = ap.parse_args()

    if args.verify:
        return run_verify(args.dirs, args.jobs)

    if args.decompress:
        todo = [os.path.join(root, name)
                for d in args.dirs for root, _, files in os.walk(d)
                for name in sorted(files) if name.endswith(".csv.zst")]
        print(f"{len(todo)} cell(s) to decompress")
        if args.dry_run:
            return 0
        failures = 0
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            for path, (size, err) in zip(todo, pool.map(decompress, todo)):
                if err:
                    failures += 1
                    print(f"  FAILED {show(path)}: {err}", file=sys.stderr)
                else:
                    print(f"  {size / 2**20:8.1f} MiB  {show(path)[:-4]}")
        return 1 if failures else 0

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

# Reproduction scripts

Scripts that drive the built binaries to regenerate the datasets behind the
paper's figures. They exist so that a result can be traced back to a command,
a commit and a seed, rather than to somebody's shell history.

Build first:

```sh
cmake --preset conan-release
cmake --build build -j 6
```

| script | drives | produces the data behind |
|---|---|---|
| `exhaustive_nonuniform.sh` | `exhaustive_search_projections` | Figures 2-6 |
| `exact_mates.sh` | `exact_mates` | Figure A.9 (post-processes the above) |
| `sample_heatmap.sh` | `count_twins_random` | Figure 7, heatmap panel |
| `sample_density.sh` | `count_twins_random` | Figure 7 line panel, Figure A.10 |
| `common.sh` | - | shared helpers; source it, do not run it |

Each script takes `--help`-free positional arguments documented in its header
comment, and is configured by environment variables.

## Conventions

**Everything is seeded.** `SEED` defaults to a fixed value, so sampled output
is a deterministic function of `(seed, sample index)` and does not depend on
the thread count. Line order within a file still varies under parallel writes,
so sort before diffing. The published sampled data predates `--seed` and cannot
be reproduced exactly; it can only be compared distributionally.

**Everything is capped.** `MAX_THREADS` defaults to 6. The drivers will
otherwise take every core on the machine.

**Every cell is guarded.** `MAX_SECONDS` and `MAX_RSS_MB` bound each individual
invocation. A cell that exceeds either is killed, its partial output deleted,
and the fact recorded - the run continues instead of taking the machine down
with it. macOS has no `timeout(1)` and no working `ulimit -v`, which is why
`run_guarded` in `common.sh` polls by hand.

**Runs are resumable.** A cell whose output already exists is skipped, so an
interrupted run can simply be restarted.

**Every run writes a manifest.** `RUN_MANIFEST.tsv` in the output directory
records the git commit, whether the tree was dirty, the host, thread count and
seed, then one row per cell giving status, wall seconds, peak RSS, output file
and the exact command. This is the provenance for anything generated here.

## Cost

Measured on an Apple M-series laptop, 6 threads. These are the numbers that
decide what is worth running; they are not from the paper.

**Exhaustive (n=6, k=3, min_k=2).** Cost is not monotone in `m`: the number of
distinct projections peaks at m=10 and then falls, while the search tree per
projection keeps growing.

| m | time | projections | peak RSS |
|---|---|---|---|
| 2-9 | 13s total | 1 -> 303 | < 0.2 GB |
| 10 | 94s | 342 | 1.2 GB |
| 11 | 39s | 304 | 1.2 GB |
| 12 | 71s | 245 | 1.9 GB |
| 13 | 124s | 159 | 4.2 GB |
| 14 | 395s | 94 | 8.7 GB |

Peak RSS roughly doubles per step above m=12, so m=15 is around 17 GB and m=16
around 34 GB. **m=2..14 is about 12 minutes in total**; beyond that memory, not
time, is the binding constraint. The default range stops at m=13.

**Heatmap sampling.** Cheap. The worst cell measured (k=4, n=16, m=16, 1000
samples) is 2.4s; k=3 cells are under 0.1s.

**Density sampling.** The expensive one. Cost climbs steeply with `m`:
n=9 m=20 is 7s, but n=9 m=40 had not finished after 10 minutes. The dense tail
of n=8 and n=9 is where the time goes, which is why `MAX_SECONDS` defaults to
300 here and the script abandons the rest of an `n` once a cell blows the
budget.

## Excluded cells

`sample_heatmap.sh` skips two cells by construction rather than failing on
them:

- **m > C(n,k)** is impossible - there are not that many distinct k-subsets.
  For k=4, n=6 that rules out m=16, since C(6,4)=15.
- **m == C(n,k)** is degenerate: the complete k-uniform hypergraph is the only
  sample, so any number of draws has an effective sample size of 1 and the
  heatmap value is identically 0.

Both are recorded in the manifest with those statuses.

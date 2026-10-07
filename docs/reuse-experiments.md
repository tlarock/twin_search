# Can previously-computed cells make expensive cells cheaper?

Two reuse ideas, both measured against real campaign output in October 2026,
both **closed negative for the samples that matter**. This records what was
measured, how, and the numbers, in enough detail to rebuild the experiment
without the branches that produced it.

The question was never "is there structure" - there is, and it is exact. It was
whether the structure removes the work that actually costs. It does not.

Experimental code lived on `experiment/nested-twin-dp` (commits `f856181`,
`5a8368c`, `9424dab`, `d44c5b1`) and may have been deleted. The reusable part,
`analysis/python/nested_twin_reuse.py`, needs no C++ changes and is kept.

---

## 0. The one lesson worth carrying forward

**Check the space cost before building the time optimisation.** Both ideas were
pursued on the strength of a promising *fraction* (how much output could be
reused, how many states collapse) without first asking what the reuse would
cost to store or whether it touched the dominant phase. In both cases a few
lines of arithmetic, done up front, would have closed the question.

---

## 1. The structural result (this part is true and exact)

`sample_uniform_random` draws hyperedges from a single
`std::uniform_int_distribution` and inserts them into a `std::set` until it
holds `m` distinct ones. The generator is seeded from `(seed, index, attempt)` -
**`m` is not in the stream** - so for a fixed index the draw sequence does not
depend on `m`, and the hyperedge SET for `m` is a prefix-set of the one for
`m+1`:

    H_{m+1}(i) = H_m(i) + {e}

(Both must cover all `n` nodes or the generator restarts the whole draw; at
campaign densities it never does. Verified: 0 violations in 500 samples for
m=34 vs m=38.)

Projections are additive, so `P_{m+1} = P_m + D(e)`, and the twin sets
decompose **exactly**:

    T(P_{m+1}) = { G' + {e} : G' in T(P_m), e not in G' }        (A)
               U { G in T(P_{m+1}) : e not in G }                (B)

(A) is a bijection: if `e` is in `G` then `P(G \ {e}) = P_m`, and conversely any
twin of `P_m` not already containing `e` extends to one of `P_{m+1}`.

**Verified on real output, not argued.** Equality of
`|{G in T(P_{m+1}) : e in G}|` and `|{G' in T(P_m) : e not in G'}|`:

| cells | samples | twins checked | result | reuse pooled | reuse median |
|---|---|---|---|---|---|
| 34->35 | 500 | 2.35M | holds | 57.9% | 83.3% |
| 35->36 | 500 | 2.73M | holds | 62.3% | 79.9% |
| 36->37 | 500 | 3.46M | holds | 57.3% | 80.5% |
| 37->38 | 463 | 2.05M | holds | 65.0% | 82.6% |
| 47->48 | 452 | 1.07M | holds | 67.3% | 85.8% |

2,415 samples, 11.7M twins, no exceptions.

### Method: offline, zero compute

No new searches are needed - the twin lists are already section 4 of every
result row.

1. `count_twins_random ... --seed <S> --dry-run | grep '^FP '` emits
   `FP i,max_log_width,num_edges,num_cliques,<projection>` for every index.
   **Run this on Linux/libstdc++** - see `docs`-adjacent note in
   `src/twin_search/count_twins_random.cpp`; `std::uniform_int_distribution` is
   implementation-defined and macOS draws different hypergraphs.
2. Rows carry no index. Recover it by rebuilding each row's projection from its
   FIRST twin (every twin has the row's projection by definition) and matching
   against the fingerprints. This is what `resume_missing.py` already does.
3. `e` is recoverable as the projection difference: adding one k=3 hyperedge
   raises exactly the three pairs it spans by one. Reject anything else as
   non-nested rather than comparing unrelated hypergraphs.
4. A twin contains `e` iff the string `|a:b:c|` occurs in `"|" + twin + "|"`.
   No parsing needed.

`analysis/python/nested_twin_reuse.py` implements all of this.

---

## 2. Why the reuse does not pay (idea 1: cross-m DP)

A DP would read (A) from disk and search only for (B). Measured with a probe
that runs each cell twice - once normally, once with `e` forbidden - and charges
the DP only for what it genuinely skips:

    dp_cost = traversal(e forbidden) + mates(full) + iso(full)

This is an upper bound on the benefit: it charges nothing for reading,
filtering and re-inflating the previous cell's twins, and the two pairwise
phases are charged in full because `G' + {e}` has a *different line graph* from
`G'`, so neither phase can reuse anything.

**Result: 1.08x aggregate over 27 real samples; 7.5% of time removed.** Worse
as samples get more expensive: 1.40x over the cheapest third, 1.05x over the
costliest.

Three separate reasons, each measured:

1. **Traversal is not where the time goes.** Pooled over that sample set:
   traversal 17.1%, line graphs + mates 10.2%, bipartites + isomorphism filter
   72.7%. The decomposition removes traversal only.
2. **Leaves are not tree nodes.** Forbidding one clique removed 43% of
   traversal where the twin-count fraction predicts 68% - an overstatement of
   1.6x. Backtracking pays for internal nodes, including dead ends, which
   forbidding a clique does not prune.
3. **A constant factor cannot outrun an exponential.** Per-sample cost growth
   is 1.16x per `m` typically, and **1.97x per m for the samples that end up
   expensive**. A 1.5x speedup therefore buys `log(1.5)/log(1.97) = 0.6` of one
   `m` step. Pulling a 100-hour sample under a 24-hour wall needs 4.2x; a
   10-day sample needs 10x.

### The population split matters, and the aggregate hides it

Pooling is misleading here. Split by what the sample spends its time on:

| population | DP speedup |
|---|---|
| pairwise-bound (large twin sets) | ~1.0x - structurally untouchable |
| traversal-bound (straggler archetype) | 1.17-3.45x, median ~1.5x |

So the DP *is* worth ~1.5x on the right population. It is still not 4x.

### Availability is a non-issue IF you commit to the chain

Initially this looked fatal: 5 of the 6 samples unfinished at m=48 were also
unfinished at m=47, so the DP's input was missing exactly where wanted. **That
is an artefact of computing cells independently, not a property of the
strategy.** Walking `m` upward always gives `j=1`. The decay table below only
applies when the chain is broken and you must skip:

| j (steps back) | pooled reuse | median/sample |
|---|---|---|
| 1 | 57.9% | 83.3% |
| 2 | 36.9% | 55.2% |
| 3 | 21.7% | 34.9% |
| 4 | 16.2% | 24.0% |

About x0.65 per step; reaching back 9-10 steps leaves ~1%. Only twins
containing ALL `j` new hyperedges are in bijection with cell `m`; partial
subsets need the intermediate cells.

Committing to the chain also **serialises the campaign** - cell `m+1` cannot
start until `m` finishes for that index - which given the QOS concurrency
limits is a large practical cost on its own.

### Method: the probe

`src/twin_search/twin_dp_probe.cpp`, plus two additions to `TwinSearch`:

* `std::set<int> forbidden_cnodes` - filtered in `get_filtered_neighbors`, which
  is the ONE place both `process_item` and `parallel_process_item` get their
  candidates, so filtering there covers both paths.
* `ms_traversal` / `ms_mates` / `ms_iso` - wall time of the three phases of
  `parallel_search`, which already has clean phase boundaries. A single total
  hides the whole point.

Per index the probe redraws `H_m` and `H_{m+1}` with the driver's own
`seed_seq{seed, i, attempt}`, confirms nesting, finds `e`, looks up its
clique-node via `fact.rev_node_map.at(e)`, then runs the full and forbidden
searches. It asserts `|T_forbidden| + |{G : e in G}| == |T_full|` every time
("split-ok"), which catches a filter applied in the wrong place.

---

## 3. Why memoising the search does not pay either (idea 2)

Stragglers are traversal-bound (section 4), so the obvious alternative is to
memoise the backtracking search itself - no previous cell needed, and it
attacks dead-end subtrees directly.

### The state is the residual projection - this was verified

The edge execution order is fixed (`default_edge_execution_order()` is just
`0..num_edge_nodes-1`), and when an edge is chosen exactly `proj_rem(e)`
cliques covering it are taken, driving it to zero and never back up. So the
next edge to satisfy is a function of `proj_rem` alone.

The non-obvious part: `get_filtered_neighbors` ALSO excludes cliques already in
the partial hypergraph, which would make the state `(proj_rem, used-set)` and
memoisation unsound. The claim is that this exclusion is redundant - a clique
was chosen to satisfy some edge, that edge is now zero, so re-selecting it
drives that entry negative and `add_to_stack` rejects it anyway.

**Tested, not assumed.** An `allow_repeat_cliques` flag drops the exclusion; the
resulting twin sets were identical on all 10 traversal-dominated samples tried,
even though the exclusion fires 0.8M-30M times per sample. The residual IS the
state.

(Incidental finding, not pursued: a check that fires tens of millions of times
per sample and never changes an answer is a linear scan over the partial
hypergraph per candidate. Whether removing it is a net win depends on whether
`add_to_stack`'s rejection is cheaper than the scan - unmeasured.)

### It still fails, on memory

Collapse ratio = nodes visited / distinct residuals = the ceiling on the time
saving. Measured across a 25x range of sample sizes:

| nodes visited | collapse |
|---|---|
| 546 K | 1.82x |
| 579 K | 3.54x |
| 738 K | 1.92x |
| 792 K | 3.24x |
| 993 K | 1.40x |
| 1.09 M | 2.74x |
| 1.71 M | 2.28x |
| 2.03 M | 2.39x |
| 2.99 M | 1.29x |
| 3.13 M | 1.90x |
| 6.30 M | 3.33x |
| 7.36 M | 2.88x |
| 7.60 M | 2.61x |
| 18.6 M | 4.92x |

Noisy, roughly flat at 2-3x, no trend with size.

A memo needs one entry per DISTINCT residual, and `distinct = visited /
collapse`. At collapse 2.5 the table is ~40% the size of the work being
avoided:

| sample | nodes visited | memo entries | memo @16 B |
|---|---|---|---|
| the 18.6M sample | 1.9e7 | 3.8e6 | 0.1 GB |
| a ~1.3-hour straggler | 3.4e10 | 1.4e10 | 218 GB |
| same, if collapse were 10x | 3.4e10 | 3.4e9 | 54 GB |
| same, if collapse were 100x | 3.4e10 | 3.4e8 | 5.4 GB |

A 24h-wall-busting sample is ~100x bigger again, and 16 B/entry is optimistic
once hash-table overhead is counted.

**So memoisation needs collapse around 100x to be memory-feasible, not 4x.**
Measured collapse is 2-5x and flat. Closed.

### Method: state instrumentation

`instrument_states` on `TwinSearch`, recording at every `process_item` /
`parallel_process_item` entry. Two things make it workable:

* **Per-thread accumulators** (`tbb::enumerable_thread_specific`), unioned at
  the end. Recording happens at every node, so a shared concurrent structure
  would serialise exactly the hot path being measured. Held behind a
  `shared_ptr` member because `std::atomic` and TBB's ETS are not copy-
  assignable and `count_twins_random` does `twins = TwinSearch(...)`.
* **HyperLogLog, not an exact hash set.** One 8-byte hash per distinct state per
  thread is more memory than the search; capping the set turns the answer into
  an unlabelled floor, which for a measurement asking "how many distinct states
  are there" is worse than no answer. 2^14 registers = 16 KB/thread regardless
  of cardinality, ~0.8% standard error. **Validated against 10 exact counts:
  error -1.62% to +1.19%.**

Residual hash is FNV-1a over the upper triangle. Collisions understate the
distinct count, which is the conservative direction here.

Instrument the PARALLEL path. Node counts are identical to sequential (same
edge order, same expansions - verified exactly on 10 samples), it is ~6.5x
faster, and `search()` still lacks the phased line-graph/bipartite build from
`90f05e1`, so the sequential path also costs roughly twice the peak memory.

---

## 4. Straggler anatomy (the most reusable result here)

There are two distinct populations of expensive sample, and conflating them
produces wrong conclusions.

**Pairwise phases scale as T^2** with a fitted constant of
**c ~= 1.0e-3 ms per twin-pair at 8 threads** (`(iso_ms + mates_ms) / T^2`,
spread 6e-4 to 1.6e-3 over the top ten by twin count). That single constant
settles which population a sample is in:

| | what costs | DP helps? | memo helps? |
|---|---|---|---|
| large twin set | O(T^2) iso + mates | no | no |
| **straggler** (huge tree, few twins) | traversal | ~1.5x | ~2.5x ceiling, infeasible memory |

**Stragglers are decisively traversal-bound.** m=48 index 314: 1,443 twins,
70,864 s recorded. Predicted iso+mates at that twin count is **2.1 seconds** -
0.04% of its cost. The conclusion is robust to any plausible calibration error.

Across a whole m=48 cell the split is roughly 2/3 of time pairwise-bound, 1/3
traversal-bound (model-based, using `c`, not measured per sample).

### The isomorphism filter is the other target, and it is a near no-op

Isomorphism classes / twins measured over whole cells: **1.000** (m=36, m=37),
0.996 (m=47), 0.999 (m=48). The filter does O(T^2) fingerprint comparisons plus
vf2 calls to remove ~0.1% of twins, and on pairwise-bound samples it is ~73% of
runtime. Replacing pairwise comparison with a canonical form (colour-refinement
/ 1-WL hash) computed once per twin would turn O(T^2) into O(T) + grouping; the
same applies to the mate test, where mates become pairs sharing a canonical
line-graph form. **Untested.** It does nothing for stragglers.

### Useful calibration

Recorded per-sample `runtime_ms` is measured under contention - 500 samples
share one TBB pool - so it is NOT a clean per-sample cost. The de-contention
factor measured against isolated 8-thread runs is **~15x**, but it is not
uniform (one sample was off by 3x). Use it for ordering and rough sizing, never
as a cost.

---

## 5. What remains open

* **Checkpointing the search** - the straggler problem is that a single sample
  cannot be SPLIT, not that 24 hours is too little compute. `resume_cell.sh`
  resumes at sample granularity, so a sample needing 30 hours never finishes on
  the short partition however often it is resubmitted. The search state is an
  explicit stack, so serialising the frontier is tractable. Gating question,
  which must be measured FIRST (see section 0): how large is the frontier?
* **Canonical-form isomorphism filtering** - section 4. Large, but only for the
  pairwise-bound population.
* **Search-order heuristics** - `compute_edge_execution_order` (most-constrained
  first) exists but is unused; the comment says it "performed similarly" to the
  default. That was presumably measured on typical samples. Stragglers are
  where variable ordering normally pays most in a CSP, and they were almost
  certainly not the test set.

"""Estimate the cost of regenerating the Figure 7 line-panel dataset.

The published samples-500 rows carry a per-sample runtime, so the cost of
re-running a cell can be estimated without running it. Two corrections are
applied:

  CALIBRATION. The published runs were made on other hardware with the older
  code. Measured against the 44 heatmap cells regenerated here, this machine
  plus the current code is ~2.2x faster per unit of work (median 2.17x, pooled
  2.21x, spread 1.36-4.37x). Pass --speedup to override.

  PARALLEL STRUCTURE. A cell's wall clock is not total_work/threads. The outer
  parallel_for_each saturates while samples remain, but the cell cannot finish
  before its single slowest sample does, and that sample is 11-68% of observed
  cell wall clock. So the estimate is

      wall ~= overhead * max(total_work / threads, slowest_sample)

  with overhead 1.07, the median ratio of measured wall clock to total_work/6
  across completed heatmap cells.

The published figures predate --seed, so the cell we would draw is NOT the cell
they drew. These are estimates of the DISTRIBUTION's cost, not of the specific
run. Heavy-tailed: a cell can land far from its estimate in either direction.

Usage:
    estimate_density_cost.py <published_density_dir> [--speedup 2.2]
                             [--threads 6] [--budget-hours 10]
"""
import os, sys, math, subprocess

# Figure 7 line panel, from figure-7-sampled-results.ipynb:
#     max_ns = {3: 9, 4: 8};  for n in range(6, max_ns[k]+1)
#     m_vals = list(range(1, math.comb(n, k)))
CURVES = {3: [6, 7, 8, 9], 4: [6, 7, 8]}
SAMPLES = 500
RUNTIME_FIELD = 4   # OLD row layout: n,m,unfiltered,filtered,runtime,width,cliques


def read_times(path):
    if path.endswith('.zst'):
        txt = subprocess.run(['zstd', '-dc', path],
                             capture_output=True, text=True).stdout
    else:
        txt = open(path).read()
    out = []
    for line in txt.splitlines():
        if not line:
            continue
        out.append(int(line.split(',')[RUNTIME_FIELD]))
    return out


def main(pub_dir, speedup=2.2, threads=6, overhead=1.07, budget_hours=10.0):
    rows, missing, partial = [], [], []
    for k, ns in CURVES.items():
        for n in ns:
            for m in range(1, math.comb(n, k)):
                stem = f'n-{n}_m-{m}_k-{k}_samples-{SAMPLES}.csv'
                p = os.path.join(pub_dir, stem + '.zst')
                if not os.path.exists(p):
                    p = os.path.join(pub_dir, stem)
                    if not os.path.exists(p):
                        missing.append((k, n, m)); continue
                t = read_times(p)
                if len(t) < SAMPLES:
                    partial.append((k, n, m, len(t)))
                    if not t:
                        continue
                    # scale up to a full 500 so the estimate is not optimistic
                    scale = SAMPLES / len(t)
                else:
                    scale = 1.0
                total = sum(t) / 1000.0 * scale
                slowest = max(t) / 1000.0
                est = overhead * max(total / speedup / threads, slowest / speedup)
                rows.append((est, k, n, m, total, slowest))
    rows.sort()

    print(f'Figure 7 line panel: {sum(len(range(1, math.comb(n,k))) for k,ns in CURVES.items() for n in ns)} cells '
          f'({len(rows)} estimable, {len(missing)} with no published data)')
    print(f'speedup {speedup}x, {threads} threads, overhead {overhead}\n')
    print(f'{"#":>4} {"k":>2}{"n":>3}{"m":>4} {"est":>10} {"cumulative":>12}   {"limited_by":>10}')
    cum = 0.0
    budget_idx = None
    for i, (est, k, n, m, total, slowest) in enumerate(rows, 1):
        cum += est
        if budget_idx is None and cum > budget_hours * 3600:
            budget_idx = i
        lim = 'tail' if slowest / speedup > total / speedup / threads else 'throughput'
        print(f'{i:>4} {k:>2}{n:>3}{m:>4} {fmt(est):>10} {fmt(cum):>12}   {lim:>10}')
    print()
    if budget_idx:
        print(f'=> {budget_idx-1} cells fit in {budget_hours:g}h; the {budget_idx}th crosses it.')
    for h in (1, 2, 4, 8, 12, 24):
        c = 0.0; n_fit = 0
        for est, *_ in rows:
            if c + est > h * 3600: break
            c += est; n_fit += 1
        print(f'   {h:>3}h -> {n_fit:>3} of {len(rows)} cells')
    if partial:
        print(f'\npublished cells with <{SAMPLES} rows (estimate scaled up): {partial}')
    if missing:
        print(f'\nno published data, cost unknown ({len(missing)}): '
              f'{[f"k{k}n{n}m{m}" for k,n,m in missing]}')

    # Per-curve view. A curve is one line on the figure; a partial curve is
    # still plottable, so what matters per curve is how much of it is cheap.
    print(f'\n{"curve":>8} {"cells":>6} {"<1m":>5} {"<10m":>6} {"<1h":>5} '
          f'{"<8h":>5} {">8h":>5} {"full curve":>12} {"cheap 90%":>11}')
    for k, ns in CURVES.items():
        for n in ns:
            c = sorted(r[0] for r in rows if r[1] == k and r[2] == n)
            if not c:
                continue
            tot = sum(c)
            cheap = sum(c[:max(1, int(len(c) * 0.9))])
            print(f'{"k="+str(k)+" n="+str(n):>8} {len(c):>6} '
                  f'{sum(1 for x in c if x<60):>5} {sum(1 for x in c if x<600):>6} '
                  f'{sum(1 for x in c if x<3600):>5} {sum(1 for x in c if x<8*3600):>5} '
                  f'{sum(1 for x in c if x>=8*3600):>5} {fmt(tot):>12} {fmt(cheap):>11}')

    # The ordered run list. m=1 and the densest cells have no published data
    # but are structurally trivial (m=1 is a single hyperedge; m=C(n,k)-1 is
    # one edge short of complete, so the twin set is forced), so they go first.
    out = os.environ.get('EMIT_LIST')
    if out:
        with open(out, 'w') as f:
            f.write('# ordered fastest-first; est_seconds is a DISTRIBUTION estimate\n')
            f.write('k\tn\tm\test_seconds\tcumulative_seconds\tnote\n')
            for k, n, m in missing:
                f.write(f'{k}\t{n}\t{m}\t-1\t-1\tno published data, expected trivial\n')
            c = 0.0
            for est, k, n, m, total, slowest in rows:
                c += est
                f.write(f'{k}\t{n}\t{m}\t{est:.1f}\t{c:.1f}\t\n')
        print(f'\nordered list written to {out}')


def fmt(s):
    if s < 90: return f'{s:.0f}s'
    if s < 5400: return f'{s/60:.1f}m'
    return f'{s/3600:.1f}h'


if __name__ == '__main__':
    a = sys.argv[1:]
    pub = a[0]
    sp = float(a[a.index('--speedup')+1]) if '--speedup' in a else 2.2
    th = int(a[a.index('--threads')+1]) if '--threads' in a else 6
    bh = float(a[a.index('--budget-hours')+1]) if '--budget-hours' in a else 10.0
    main(pub, sp, th, 1.07, bh)

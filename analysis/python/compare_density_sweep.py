#!/usr/bin/env python3
"""Compare the re-run density sweep against an older sampled run, statistically.

Companion to gram_mates' compare_sampled_figure7.py, which hardcodes 1000
samples in BOTH arms. The density cells do not match that: the older run used
100 samples per cell and the re-run uses 500, so the arms have different n and
the filenames differ. Everything else - the statistic, the traps - is the same.

Headline statistic, as in Figure 7: p = #{samples with total_filtered > 1} / n.

The published samples predate --seed, so this is necessarily DISTRIBUTIONAL:
two independent draws from the same generative model, one per code version.

The two traps from the companion script apply unchanged:
  1. Do NOT test the per-cell p-values for uniformity. Fisher's exact is
     discrete, so under the null its p-values are stochastically LARGER than
     U(0,1). Count p < 0.05 against its expectation instead.
  2. Do NOT compare max_width multisets for equality. Different random
     hypergraphs have different widths; that says nothing about the code.
"""
import os, re, sys
import numpy as np
from scipy import stats

PAT = re.compile(r'^n-(\d+)_m-(\d+)_k-(\d+)_samples-(\d+)\.csv$')


def parse(path):
    mw, tw, tf, tm = [], [], [], []
    for line in open(path):
        h, t = line.find('/'), line.rfind('/') + 1
        f = line[:h].split(',')
        mw.append(int(f[3])); tm.append(int(f[6]))
        mid = line[h + 1:t - 1].split('/')
        unf = filt = 0
        for i in range(0, len(mid) - 1, 2):
            a, b = mid[i + 1].split(',')
            unf += int(a); filt += int(b)
        tw.append(unf); tf.append(filt)
    return dict(max_width=np.array(mw), total_twins=np.array(tw),
                total_filtered=np.array(tf), total_mates=np.array(tm))


def bh(p):
    p = np.asarray(p, float); order = np.argsort(p); n = len(p)
    adj = np.empty(n); prev = 1.0
    for rank, idx in enumerate(reversed(order), 1):
        prev = min(prev, p[idx] * n / (n - rank + 1)); adj[idx] = prev
    return adj


def index(d):
    """-> {(k,n,m): (path, samples)}, keeping the LARGEST sample count.

    A cell can exist at several sample counts in the same directory (4 do in
    results/increasing_density). Taking whichever os.listdir yielded last made
    the comparison depend on directory order.
    """
    out = {}
    for f in sorted(os.listdir(d)):
        mo = PAT.match(f)
        if mo:
            n, m, k, s = map(int, mo.groups())
            if (k, n, m) not in out or s > out[(k, n, m)][1]:
                out[(k, n, m)] = (os.path.join(d, f), s)
    return out


def main(old_dir, new_dir, k_want=None):
    old, new = index(old_dir), index(new_dir)
    common = sorted(set(old) & set(new))
    if k_want is not None:
        common = [c for c in common if c[0] == k_want]
    rows, pooled = [], dict(om=0, on=0, nm=0, nn=0)
    for key in common:
        (op, os_), (np_, ns) = old[key], new[key]
        o, w = parse(op), parse(np_)
        if len(o['total_filtered']) != os_ or len(w['total_filtered']) != ns:
            print(f"  SKIP {key}: row count != filename", file=sys.stderr); continue
        om, nm = int((o['total_filtered'] > 1).sum()), int((w['total_filtered'] > 1).sum())
        on, nn = len(o['total_filtered']), len(w['total_filtered'])
        for a, b in (('om', om), ('on', on), ('nm', nm), ('nn', nn)): pooled[a] += b
        fisher = stats.fisher_exact([[om, on - om], [nm, nn - nm]])[1]
        mwu = stats.mannwhitneyu(o['total_twins'], w['total_twins']).pvalue
        ks = stats.ks_2samp(o['max_width'], w['max_width']).pvalue
        rows.append((key, om / on, nm / nn, on, nn, fisher, mwu, ks))

    print(f"{'k':>2} {'n':>2} {'m':>3} {'p_old':>7} {'p_new':>7} {'delta':>8} "
          f"{'n_old':>6} {'n_new':>6} {'fisher':>8} {'MWU':>8} {'KS':>8}")
    for (k, n, m), po, pn, on, nn, f, u, s in rows:
        print(f"{k:>2} {n:>2} {m:>3} {po:>7.3f} {pn:>7.3f} {pn-po:>+8.3f} "
              f"{on:>6} {nn:>6} {f:>8.4f} {u:>8.4f} {s:>8.4f}")
    if not rows:
        print("no overlapping cells"); return
    fp = [r[5] for r in rows]
    print(f"\ncells {len(rows)}   max |delta p| {max(abs(r[2]-r[1]) for r in rows):.4f}"
          f"   mean delta {np.mean([r[2]-r[1] for r in rows]):+.4f}")
    print(f"Fisher p<0.05: {sum(p < .05 for p in fp)}  (expected ~{0.05*len(rows):.1f} by chance)"
          f"; after BH: {int((bh(fp) < .05).sum())}")
    print(f"Mann-Whitney on total_twins p<0.05: {sum(r[6] < .05 for r in rows)}")
    print(f"KS on max_width p<0.05: {sum(r[7] < .05 for r in rows)}")
    po, pn = pooled['om']/pooled['on'], pooled['nm']/pooled['nn']
    print(f"\nNAIVE pooling ({pooled['on']} old / {pooled['nn']} new samples)")
    print(f"  P(multi-iso) old {po:.4f}   new {pn:.4f}   difference {pn-po:+.4f}")
    print("  *** CONFOUNDED - do not quote. p varies from ~0 to ~1 across cells and")
    print("      the arms weight those cells differently (100 vs 500 samples each,")
    print("      and not the same set of cells at each count). This is Simpson's")
    print("      paradox territory: the companion script could pool safely only")
    print("      because both of its arms had exactly 1000 samples in every cell.")

    # Cochran-Mantel-Haenszel: tests a common odds ratio ACROSS strata, so each
    # cell is compared only against itself. This is the pooled test that means
    # something here.
    num = den = 0.0; a_sum = e_sum = v_sum = 0.0
    for (k, n_, m_), po_c, pn_c, on, nn, *_ in rows:
        a = round(po_c*on); b = on-a; c = round(pn_c*nn); d = nn-c
        N = a+b+c+d
        num += a*d/N; den += b*c/N
        a_sum += a; e_sum += (a+b)*(a+c)/N
        v_sum += (a+b)*(c+d)*(a+c)*(b+d)/(N*N*(N-1))
    orr = num/den if den else float('nan')
    chi2 = (abs(a_sum-e_sum)-0.5)**2/v_sum if v_sum else float('nan')
    pcmh = stats.chi2.sf(chi2, 1)
    print(f"\nCochran-Mantel-Haenszel, stratified by cell ({len(rows)} strata)")
    print(f"  common odds ratio (old:new) {orr:.4f}   chi2 {chi2:.2f}   p {pcmh:.3f}")
    print(f"  mean per-cell delta p {np.mean([r[2]-r[1] for r in rows]):+.4f}")
    print("\nA null result BOUNDS the discrepancy, it does not prove zero.")


if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0], a[1], int(a[2]) if len(a) > 2 else None)

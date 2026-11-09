# Self-contained demo: plant trend-preserving biclusters, add noise, compare strict LCS with approximate LCS.
import os
import time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from bicluster import find_biclusters, match_score


def plant(n_rows=240, n_cols=40, k=3, rows_per=25, cols_per=12, noise=0.1, seed=0):
    """Background N(0,1). Each bicluster: its rows share one column ordering (different scales/offsets)."""
    rng = np.random.default_rng(seed)
    M = rng.normal(0, 1, (n_rows, n_cols))
    truth = []
    rperm, cperm = rng.permutation(n_rows), rng.permutation(n_cols)
    for b in range(k):
        rows = rperm[b * rows_per:(b + 1) * rows_per]
        cols = rng.choice(n_cols, cols_per, replace=False)  # biclusters may overlap in columns
        base = np.sort(rng.normal(0, 1, cols_per))
        for r in rows:
            M[r, cols] = rng.uniform(0.5, 2) * base + rng.normal(0, 0.5) + rng.normal(0, noise, cols_per)
        truth.append((rows, cols))
    return M, truth


def main():
    os.makedirs("figures", exist_ok=True)
    os.makedirs("results", exist_ok=True)
    rows = []
    for noise in [0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4]:
        for rep in range(4):
            M, truth = plant(noise=noise, seed=rep)
            for name, tol in [("strict LCS (tol=0)", 0.0), ("approximate LCS (tol=0.1)", 0.1)]:
                t = time.time()
                found = find_biclusters(M, tol=tol, seed=rep)
                rows.append({"noise_sd": noise, "rep": rep, "method": name,
                             "recovery": match_score(truth, found), "relevance": match_score(found, truth),
                             "n_found": len(found), "seconds": time.time() - t})
    raw = pd.DataFrame(rows)
    table = raw.groupby(["noise_sd", "method"])[["recovery", "relevance", "n_found"]].mean().round(3).reset_index()
    table.to_csv("results/summary.csv", index=False)

    M, truth = plant(noise=0.2, seed=0)
    found = find_biclusters(M, tol=0.1)
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.2), gridspec_kw={"width_ratios": [1, 1, 1.2]})
    r0, c0 = truth[0]
    rest_c = np.setdiff1d(np.arange(M.shape[1]), c0)
    axes[0].plot(M[np.ix_(r0[:12], np.r_[c0[np.argsort(M[r0[0], c0])], rest_c])].T, lw=1)
    axes[0].axvline(len(c0) - 0.5, color="k", ls="--")
    axes[0].set_title("A trend-preserving bicluster (noise sd 0.2)\nsame rise across its 12 columns, then noise")
    axes[0].set_xlabel("columns: bicluster columns in shared order | the rest")
    axes[0].set_ylabel("expression")

    best = max(found, key=lambda f: len(set(f[0]) & set(r0)))
    rr = np.r_[best[0], np.setdiff1d(np.arange(M.shape[0]), best[0])[:40]]
    cc = np.r_[best[1], np.setdiff1d(np.arange(M.shape[1]), best[1])]
    from bicluster import minmax_rows
    axes[1].imshow(minmax_rows(M)[np.ix_(rr, cc)], aspect="auto", cmap="viridis")
    axes[1].axhline(len(best[0]) - 0.5, color="w"); axes[1].axvline(len(best[1]) - 0.5, color="w")
    axes[1].set_title(f"Recovered: {len(best[0])} rows x {len(best[1])} cols (top-left)\nrows min-max scaled")
    axes[1].set_xlabel("columns"); axes[1].set_ylabel("rows")

    for name, c in [("strict LCS (tol=0)", "#c44e52"), ("approximate LCS (tol=0.1)", "#4c72b0")]:
        sub = table[table.method == name]
        axes[2].plot(sub.noise_sd, sub.recovery, "o-", color=c, label=f"{name}: recovery")
        axes[2].plot(sub.noise_sd, sub.relevance, "s--", color=c, alpha=0.6, label=f"{name}: relevance")
    axes[2].set_xlabel("noise sd added inside biclusters")
    axes[2].set_ylabel("match score (Jaccard, Prelic et al.)")
    axes[2].set_ylim(0, 1.02)
    axes[2].set_title("Tolerating near-ties keeps biclusters intact as noise grows")
    axes[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig("figures/demo.png", dpi=120)
    print(table.to_string(index=False))
    print(f"mean seconds per run: {raw.seconds.mean():.1f}")


if __name__ == "__main__":
    main()

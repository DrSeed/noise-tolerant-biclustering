# Trend-preserving biclustering with an approximate longest common subsequence, after NoiBic (simplified).
from bisect import bisect_right
import numpy as np


def minmax_rows(M):
    lo, hi = M.min(axis=1, keepdims=True), M.max(axis=1, keepdims=True)
    return (M - lo) / (hi - lo + 1e-12)


def _longest_chain(order_vals, cols):
    """Longest non-decreasing subsequence of order_vals; returns the matching columns."""
    tails, tail_idx, parent = [], [], [-1] * len(order_vals)
    for i, v in enumerate(order_vals):
        j = bisect_right(tails, v)
        if j == len(tails):
            tails.append(v); tail_idx.append(i)
        else:
            tails[j] = v; tail_idx[j] = i
        parent[i] = tail_idx[j - 1] if j > 0 else -1
    out, i = [], tail_idx[-1] if tail_idx else -1
    while i >= 0:
        out.append(cols[i]); i = parent[i]
    return out[::-1]


def lacs(xu, xv, cols, tol):
    """Columns on which rows u and v rise together. Values closer than `tol` count as ties, so
    noise-induced swaps of near-equal values don't break the chain (tol=0: classic LCS of orders)."""
    cols = np.asarray(cols)
    if tol > 0:
        qu, qv = np.floor(xu[cols] / tol), np.floor(xv[cols] / tol)
    else:
        qu, qv = xu[cols], xv[cols]
    order = np.lexsort((qv, qu))
    return _longest_chain(list(qv[order]), list(cols[order]))


def find_biclusters(M, n_biclusters=3, tol=0.1, eta=0.85, n_groups=4, min_rows=5, seed=0):
    """Seeds = row pairs with the longest shared trend; grow rows while pruning columns to the
    common trend; then noise-aware row and column expansion. Returns [(rows, cols), ...]."""
    rng = np.random.default_rng(seed)
    X = minmax_rows(M)
    n, m = X.shape
    allc = np.arange(m)
    seeds = []
    for group in np.array_split(rng.permutation(n), n_groups):  # pairs within random row groups, as in NoiBic
        for a in range(len(group)):
            for b in range(a + 1, len(group)):
                u, v = group[a], group[b]
                seeds.append((len(lacs(X[u], X[v], allc, tol)), u, v))
    seeds.sort(reverse=True)

    found, used = [], set()
    for _, u, v in seeds:
        if len(found) == n_biclusters:
            break
        if u in used or v in used:
            continue
        rows, cols = [u, v], lacs(X[u], X[v], allc, tol)
        others = sorted((r for r in range(n) if r not in rows),
                        key=lambda r: -len(lacs(X[u], X[r], cols, tol)))
        for r in others:  # add rows, pruning columns to the trend shared by every row
            new_cols = lacs(X[u], X[r], cols, tol)
            new_cols = [c for c in cols if c in set(new_cols)]
            if min(len(rows) + 1, len(new_cols)) < min(len(rows), len(cols)):
                break
            rows.append(r); cols = new_cols
        cols = _expand_columns(X, rows, cols, allc, tol, eta)
        rows = [r for r in range(n) if r in rows or (r not in used and
                len(lacs(X[u], X[r], cols, tol)) >= eta * len(cols))]
        if len(rows) >= min_rows:
            found.append((np.array(rows), np.array(cols)))
            used.update(rows)
    return found


def _expand_columns(X, rows, cols, allc, tol, eta):
    """Insert a new column if, in at least eta of the rows, its value lands at the same place
    (+-1 position) in the bicluster's shared order."""
    cols = list(cols)
    for c in allc:
        if c in cols:
            continue
        sub = X[np.ix_(rows, cols)]
        pos = (sub < X[rows, c][:, None] - tol / 2).sum(axis=1)
        if np.mean(np.abs(pos - np.median(pos)) <= 1) >= eta:
            cols.insert(int(np.median(pos)), c)
    return cols


def match_score(A, B):
    """Prelic et al. (2006): mean over biclusters in A of the best Jaccard overlap (cells) with B."""
    if not A:
        return 0.0
    def cells(bc):
        return {(r, c) for r in bc[0] for c in bc[1]}
    sb = [cells(b) for b in B]
    return float(np.mean([max(len(ca & s) / len(ca | s) for s in sb) if sb else 0 for ca in map(cells, A)]))

"""Association tests: EMMAX-style linear mixed model, hypergeometric baseline, BH.

Model choice and defaults come from the calibration grid in docs/calibration/README.md."""

import numpy as np
import pandas as pd
from scipy import optimize, stats


def kinship(distance: pd.DataFrame) -> pd.DataFrame:
    """Genetic relatedness from a Hamming-type distance (already a squared Euclidean distance):
    double-centred -D/2 (classical MDS), mean diagonal 1. Do not square D: it breaks positive semi-definiteness."""
    D = distance.to_numpy(float)
    J = np.eye(len(D)) - 1 / len(D)
    B = -0.5 * J @ D @ J
    return pd.DataFrame(B / np.diag(B).mean(), index=distance.index, columns=distance.columns)


def pcs(K: pd.DataFrame, k: int) -> pd.DataFrame:
    """Top-k principal coordinates of a kinship matrix."""
    s, U = np.linalg.eigh(K.to_numpy())
    top = np.argsort(s)[::-1][:k]
    return pd.DataFrame(U[:, top] * np.sqrt(np.clip(s[top], 0, None)), index=K.index,
                        columns=[f"PC{i + 1}" for i in range(k)])


def inverse_normal(y: np.ndarray) -> np.ndarray:
    """Rank-based inverse normal transform; keeps the t-test calibrated when carriers are few and y is skewed."""
    return stats.norm.ppf((stats.rankdata(y) - 0.5) / len(y))


def lmm_null(K: np.ndarray, y: np.ndarray, X: np.ndarray) -> np.ndarray:
    """Fit y = X b + u + e, u ~ N(0, sg2*K), e ~ N(0, se2*I) once under the null (EMMAX, Kang 2010).
    Returns W such that W @ M whitens any strains-x-k matrix M; then OLS on W-rotated data is the LMM test."""
    s, U = np.linalg.eigh(K)
    s = np.clip(s, 0, None)
    yt, Xt = U.T @ y, U.T @ X

    def nll(log_delta):  # ML profile likelihood of delta = se2/sg2
        v = s + np.exp(log_delta)
        Xw = Xt / v[:, None]
        r = yt - Xt @ np.linalg.solve(Xt.T @ Xw, Xw.T @ yt)
        return 0.5 * (len(y) * np.log((r * r / v).sum()) + np.log(v).sum())

    delta = np.exp(optimize.minimize_scalar(nll, bounds=(-10, 10), method="bounded").x)
    return U.T / np.sqrt(s + delta)[:, None]


def assoc(G: np.ndarray, y: np.ndarray, X: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """OLS y ~ X + g for each column g of G, all at once (Frisch-Waugh-Lovell).
    X must include the intercept. Returns beta, se, p; NaN where g is constant after adjustment."""
    Q, _ = np.linalg.qr(X)
    yr = y - Q @ (Q.T @ y)
    Gr = G - Q @ (Q.T @ G)
    ss = (Gr**2).sum(axis=0)
    df = len(y) - X.shape[1] - 1
    with np.errstate(divide="ignore", invalid="ignore"):
        ss = np.where(ss > 1e-12, ss, np.nan)
        beta = Gr.T @ yr / ss
        se = np.sqrt(((yr**2).sum() - beta**2 * ss) / df / ss)
        p = 2 * stats.t.sf(np.abs(beta / se), df)
    return beta, se, p


def hypergeom_p(carrier: np.ndarray, y: np.ndarray, sd: float = 1.5) -> float:
    """Original Mutation Hunter test: are carriers enriched among strains with y > mean + sd*SD?"""
    extreme = y > y.mean() + sd * y.std(ddof=1)
    return float(stats.hypergeom.sf((carrier & extreme).sum() - 1, len(y), extreme.sum(), carrier.sum()))


def bh(p) -> np.ndarray:
    """Benjamini-Hochberg q-values; NaN p-values stay NaN and do not count as tests."""
    p = np.asarray(p, float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if ok.any():
        q[ok] = stats.false_discovery_control(p[ok])
    return q

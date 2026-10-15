import numpy as np
import pandas as pd
from scipy import stats as sps

from mutation_hunter import stats


def test_assoc_matches_full_ols():
    rng = np.random.default_rng(1)
    n = 60
    C = rng.normal(size=(n, 2))
    g = (rng.random(n) < 0.2).astype(float)
    y = 0.5 * g + C @ [0.3, -0.2] + rng.normal(size=n)
    beta, se, p = stats.assoc(g[:, None], y, np.column_stack([np.ones(n), C]))
    X = np.column_stack([np.ones(n), C, g])
    coef, rss, *_ = np.linalg.lstsq(X, y, rcond=None)
    se_full = np.sqrt(rss[0] / (n - 4) * np.linalg.inv(X.T @ X)[3, 3])
    assert np.isclose(beta[0], coef[3]) and np.isclose(se[0], se_full)
    assert np.isclose(p[0], 2 * sps.t.sf(abs(coef[3] / se_full), n - 4))


def test_assoc_constant_column_is_nan():
    beta, se, p = stats.assoc(np.zeros((10, 1)), np.arange(10.0), np.ones((10, 1)))
    assert np.isnan(p[0])


def test_hypergeom_matches_original_r_script():
    y = np.array([0, 0, 0, 0, 0, 0, 0, 0, 10, 10.0])  # 2 extreme strains
    carrier = np.array([0, 0, 0, 0, 0, 0, 0, 1, 1, 1], bool)
    # R: phyper(k-1, N, M-N, n, lower.tail=FALSE) with M=10, N=2, n=3, k=2
    assert np.isclose(stats.hypergeom_p(carrier, y), sps.hypergeom.sf(1, 10, 2, 3))


def test_bh_skips_nan():
    q = stats.bh([0.01, np.nan, 0.04])
    assert np.isnan(q[1]) and np.allclose(q[[0, 2]], [0.02, 0.04])


def test_lmm_with_identity_kinship_equals_ols():
    rng = np.random.default_rng(2)
    n = 50
    g = (rng.random(n) < 0.3).astype(float)[:, None]
    y = g[:, 0] + rng.normal(size=n)
    ones = np.ones((n, 1))
    W = stats.lmm_null(np.eye(n), y, ones)  # K = I: no structure, W only rescales
    assert np.isclose(stats.assoc(W @ g, W @ y, W @ ones)[2][0], stats.assoc(g, y, ones)[2][0])


def test_kinship_relatives_share_more():
    d = pd.DataFrame([[0, .1, 1, 1], [.1, 0, 1, 1], [1, 1, 0, .1], [1, 1, .1, 0]], dtype=float)
    K = stats.kinship(d).to_numpy()
    assert K[0, 1] > 0 > K[0, 2] and np.isclose(np.diag(K).mean(), 1)


def test_inverse_normal_is_rank_based_and_symmetric():
    z = stats.inverse_normal(np.array([0.1, 5.0, 0.2, 100.0]))
    assert list(np.argsort(z)) == [0, 2, 1, 3] and np.isclose(z.sum(), 0)


def test_kinship_from_hamming_distance_is_positive_semidefinite():
    rng = np.random.default_rng(3)
    X = rng.integers(0, 2, (30, 200))
    d = pd.DataFrame((X[:, None, :] != X[None, :, :]).mean(axis=2))  # Hamming fraction, like the public matrices
    assert np.linalg.eigvalsh(stats.kinship(d).to_numpy()).min() > -1e-9

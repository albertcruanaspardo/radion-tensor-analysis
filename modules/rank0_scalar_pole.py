import numpy as np
from scipy import stats

def evaluate_rank0(m_true=602.50, gamma_true=1.739, n_obs=25, b_exp=2.59):
    p_val = float(stats.poisson.sf(n_obs - 1, b_exp))
    z_score = float(stats.norm.ppf(1.0 - p_val)) if p_val > 0 else 8.21
    delta_chi2 = float(z_score**2)
    return {
        "rank": 0,
        "name": "Conformal Scalar Pole (ZZ -> 4l)",
        "delta_chi2": delta_chi2,
        "z_score": z_score,
        "p_val": p_val
    }

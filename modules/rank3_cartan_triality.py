import numpy as np
from scipy import stats

def evaluate_rank3(a3_measured=0.1873, a3_err=0.0637):
    z_score = float(a3_measured / a3_err)
    delta_chi2 = float(z_score**2)
    p_val = float(stats.norm.sf(abs(z_score)) * 2.0)
    return {
        "rank": 3,
        "name": "Cartan SO(8) Chiral Torsion",
        "delta_chi2": delta_chi2,
        "z_score": z_score,
        "p_val": p_val
    }

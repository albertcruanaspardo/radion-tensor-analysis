import numpy as np
from scipy import stats, optimize

def evaluate_rank1():
    datasets = {
        "Cryo_Cavities": {"obs": 0.40e-18, "err": 1.20e-18, "sm": 0.0, "scale": 1e18},
        "LHC_A_FB": {"obs": 0.0612, "err": 0.0043, "sm": 0.0618, "scale": 1.0},
        "CMS_CM_Triad": {"obs": 0.00518, "err": 0.07254, "sm": 0.0, "scale": 1.0},
        "IceCube_Zenith": {"obs": 0.00140, "err": 0.01180, "sm": 0.0, "scale": 1.0}
    }
    y_obs = np.array([d["obs"] * d["scale"] for d in datasets.values()])
    y_sm = np.array([d["sm"] * d["scale"] for d in datasets.values()])
    sigmas = np.array([d["err"] * d["scale"] for d in datasets.values()])

    chi2_sm = float(np.sum(((y_obs - y_sm) / sigmas)**2))
    p_val_sm = float(stats.chi2.sf(chi2_sm, df=len(y_obs)))

    sens = np.array([1.0, 0.45, 0.85, 0.30])
    def chi2_func(eta):
        return np.sum(((y_obs - (y_sm + sens * eta)) / sigmas)**2)

    res = optimize.minimize_scalar(chi2_func, bounds=(-1.0, 1.0), method="bounded")
    chi2_min = float(res.fun)
    delta_chi2 = max(0.0, chi2_sm - chi2_min)
    z_score = float(np.sqrt(delta_chi2))

    return {
        "rank": 1,
        "name": "Foliation Confinement (Vector Nullity)",
        "delta_chi2": delta_chi2,
        "z_score": z_score,
        "p_val": p_val_sm
    }

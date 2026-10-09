import numpy as np
from scipy import stats, optimize

def evaluate_rank4():
    channels = {
        "CMS_VBS_ssWW": {"obs": 47, "bkg": 36.8, "lin": 1.45, "quad": 0.082, "sys": 0.12},
        "ATLAS_VBS_ssWW": {"obs": 52, "bkg": 41.2, "lin": 1.55, "quad": 0.088, "sys": 0.11},
        "CMS_VBS_WZ": {"obs": 23, "bkg": 17.5, "lin": 0.95, "quad": 0.052, "sys": 0.14},
        "ATLAS_VBS_WZ": {"obs": 28, "bkg": 21.0, "lin": 1.05, "quad": 0.058, "sys": 0.13},
        "LHC_ZZ_4l_VBS": {"obs": 8, "bkg": 4.9, "lin": 0.62, "quad": 0.035, "sys": 0.09},
        "IceCube_193TeV": {"obs": 114, "bkg": 110.2, "lin": 0.85, "quad": 0.015, "sys": 0.15}
    }
    def nll(c):
        val = 0.0
        for ch in channels.values():
            mu = max(ch["bkg"] + ch["lin"] * c + ch["quad"] * (c**2), 1e-6)
            var = mu + (ch["sys"] * ch["bkg"])**2
            val += 0.5 * ((ch["obs"] - mu)**2 / var)
        return val

    nll_sm = nll(0.0)
    res = optimize.minimize_scalar(nll, bounds=(-2.0, 15.0), method="bounded")
    delta_chi2 = float(2.0 * (nll_sm - res.fun))
    z_score = float(np.sqrt(max(0.0, delta_chi2)))
    p_val = float(stats.chi2.sf(delta_chi2, df=1))
    return {
        "rank": 4,
        "name": "Dimension-8 SMEFT & VBS Unitarization",
        "delta_chi2": delta_chi2,
        "z_score": z_score,
        "p_val": p_val
    }

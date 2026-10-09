import numpy as np
from scipy import stats, optimize

def evaluate_rank2(v=246.22, m_h=125.25, m_phi=602.50, lambda_r=3000.0):
    gamma = v / lambda_r
    kappas = [
        (1.012, 0.038, "V"), (0.993, 0.040, "V"),
        (0.970, 0.060, "F"), (0.930, 0.070, "F"),
        (0.985, 0.055, "F"), (1.020, 0.048, "L")
    ]
    hh_channels = [
        (1.45, 0.38, "mu"), (1.30, 0.45, "mu"), (1.15, 0.55, "mu"),
        (2.85, 1.10, "sig_bbyy"), (14.2, 5.80, "sig_bbta")
    ]
    def predict(xi):
        num = 12.0 * xi * gamma * (1.0 + (m_h / m_phi)**2)
        den = 1.0 - (m_h / m_phi)**2 - 6.0 * (xi**2) * (gamma**2) * (1.0 - (m_h / m_phi)**2)
        theta = 0.5 * np.arctan(num / max(den, 1e-6))
        cos_th, sin_th = np.cos(theta), np.sin(theta)
        kV = cos_th - gamma * sin_th
        kF = cos_th + gamma * sin_th
        kL = 1.28 * kV - 0.28 * kF
        c_curv = 1.0 + 2.0 * xi * (1.0 + 2.0 * (m_h / m_phi)**2)
        mu_hh = (kV**2) * (1.0 + 0.72 * xi + 0.15 * xi**2)
        sig_bbyy = 2.10 * (c_curv**2) * (cos_th**2)
        sig_bbta = 11.5 * (c_curv**2) * (cos_th**2)
        return [kV, kV, kF, kF, kF, kL], [mu_hh, mu_hh, mu_hh, sig_bbyy, sig_bbta]

    def chi2_global(params):
        xi = params[0]
        pk, phh = predict(xi)
        c2_k = sum(((obs - p) / err)**2 for (obs, err, _), p in zip(kappas, pk))
        c2_hh = sum(((obs - p) / err)**2 for (obs, err, _), p in zip(hh_channels, phh))
        return c2_k + c2_hh

    c2_sm = sum(((obs - 1.0) / err)**2 for obs, err, _ in kappas) +             sum(((obs - 1.0) / err)**2 for obs, err, t in hh_channels if t == "mu") +             sum(((obs - 0.0) / err)**2 for obs, err, t in hh_channels if t.startswith("sig"))

    res = optimize.minimize(chi2_global, x0=[0.06], bounds=[(-0.4, 0.4)], method="L-BFGS-B")
    delta_chi2 = float(c2_sm - res.fun)
    z_score = float(np.sqrt(max(0.0, delta_chi2)))
    p_val = float(stats.chi2.sf(delta_chi2, df=1))
    return {
        "rank": 2,
        "name": "Ricci Curvature Mixing & Di-Higgs",
        "delta_chi2": delta_chi2,
        "z_score": z_score,
        "p_val": p_val
    }

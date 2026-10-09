#!/usr/bin/env python3
"""
Octet Spectral Miner v2: Incorporating Parton Luminosity Weighting
Author: Albert Cruañas Pardo (ORCID: 0009-0007-0846-2321)
"""

import numpy as np
from scipy import optimize, stats

octeto_data = [
    {"n": 1, "lambda": 3000.0, "dm2": +84.43, "canal": "Drell-Yan / gamma-gamma", "obs": 35, "bkg": 24.65, "sys": 0.08},
    {"n": 2, "lambda": 2566.0, "dm2": +49.96, "canal": "Dijet m_jj in [2.4, 2.7] TeV", "obs": 182, "bkg": 158.30, "sys": 0.06},
    {"n": 3, "lambda": 2194.8, "dm2": +57.37, "canal": "c-tagged dijet / Z+c", "obs": 78, "bkg": 63.80, "sys": 0.09},
    {"n": 4, "lambda": 1877.4, "dm2": +31.29, "canal": "Tau+ Tau- (hadronic boosted)", "obs": 44, "bkg": 34.20, "sys": 0.10},
    {"n": 5, "lambda": 1605.8, "dm2": +15.47, "canal": "b-bbar boosted pair", "obs": 310, "bkg": 284.10, "sys": 0.05},
    {"n": 6, "lambda": 1373.5, "dm2":  +6.01, "canal": "Monojet + MET (pT > 1.2 TeV)", "obs": 520, "bkg": 496.00, "sys": 0.04},
    {"n": 7, "lambda": 1174.8, "dm2":  +0.75, "canal": "t-tbar fat-jet resonance", "obs": 415, "bkg": 402.50, "sys": 0.04},
    {"n": 8, "lambda": 1004.9, "dm2":  -3.86, "canal": "Di-Higgs hh -> 4b/bbgg", "obs": 62, "bkg": 54.10, "sys": 0.07}
]

sqrt_s = 13000.0 # GeV

# Pesado físico de señal por luminosidad diferencial partónica efectiva
s_phys = []
for row in octeto_data:
    m = row["lambda"]
    # Aproximación cinemática de luminosidad diferencial dL/dm
    tau = (m / sqrt_s)**2
    pdf_factor = ((1.0 - m / sqrt_s)**5.5) / (m / 1000.0)
    # Acoplamiento proporcional a la curvatura neta
    s_val = (abs(row["dm2"]) / 10.0) * pdf_factor * 12.5
    s_phys.append(s_val)

def nll_v2(mu):
    chi2 = 0.0
    for i, row in enumerate(octeto_data):
        b = row["bkg"]
        signo = 1.0 if row["dm2"] >= 0 else -0.3
        s = max(0.0, mu) * s_phys[i] * signo
        mu_tot = max(1e-5, b + s)
        sigma2 = mu_tot + (row["sys"] * b)**2
        chi2 += ((row["obs"] - mu_tot)**2) / sigma2
    return 0.5 * chi2

res_0 = nll_v2(0.0)
res = optimize.minimize_scalar(nll_v2, bounds=(0.01, 10.0), method="bounded")
mu_hat = float(res.x)
delta_chi2 = 2.0 * (res_0 - res.fun)
z_score = float(np.sqrt(max(0.0, delta_chi2)))
p_val = float(stats.chi2.sf(delta_chi2, df=1))

print("=" * 95)
print("MINERÍA ESPECTRAL RECALIBRADA CON LUMINOSIDAD PARTÓNICA:")
print("=" * 95)
print(f"{'Rango':<6} | {'Escala':<12} | {'Obs':<5} | {'Bkg':<7} | {'Exceso':<7} | {'S (pred Fís.)'}")
print("-" * 95)
for i, row in enumerate(octeto_data):
    exc = row["obs"] - row["bkg"]
    print(f"n = {row['n']:<2} | {row['lambda']:6.1f} GeV  | {row['obs']:<5} | {row['bkg']:7.2f} | {exc:+7.2f} | {s_phys[i]*mu_hat:6.2f}")
print("-" * 95)
print(f"[*] Modulador Global Ajustado (mu_hat): {mu_hat:.3f}")
print(f"[*] Delta Chi2 frente al SM puro:       {delta_chi2:.2f}")
print(f"[*] SIGNIFICANCIA GLOBAL DEL OCTETO:    Z = {z_score:.2f} sigma")
print(f"[*] p-value:                            p = {p_val:.3e}")
print("=" * 95)

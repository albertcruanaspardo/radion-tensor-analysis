#!/usr/bin/env python3
"""
Octet Spectral Pattern Miner for ATLAS and CMS Invariant Mass Distributions
Evaluates the correlated 8-level imprint (n = 1..8) across LHC Run 2 & Run 3 bins.
Author: Albert Cruañas Pardo (ORCID: 0009-0007-0846-2321)
"""

import numpy as np
from scipy import optimize, stats

print("=" * 95)
print("MINERÍA ESPECTRAL DEL DOBLE OCTETO (n = 1 .. 8) EN DATOS DEL LHC (ATLAS + CMS)")
print("=" * 95)

# Matriz del Octeto calculada
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

# Normalización del vector de señal teórica predicha según el balance de masa
# S_n nominal proporcional a la densidad de modulación cuántica (dm2 / lambda)
s_nom = []
for row in octeto_data:
    # Señal normalizada al orden de magnitud esperado por luminosidad L = 140 fb^-1
    s_teorica = abs(row["dm2"]) / (row["lambda"] * 0.012)
    s_nom.append(s_teorica)

# Verosimilitud de Poisson con sistemáticas gaussianas
def nll_octeto(mu):
    chi2 = 0.0
    for i, row in enumerate(octeto_data):
        b = row["bkg"]
        # El nivel 8 tiene interferencia destructiva (signo negativo en dm2)
        signo = 1.0 if row["dm2"] >= 0 else -0.4
        s = max(0.0, mu) * s_nom[i] * signo
        mu_tot = max(1e-5, b + s)
        sigma2 = mu_tot + (row["sys"] * b)**2
        chi2 += ((row["obs"] - mu_tot)**2) / sigma2
    return 0.5 * chi2

# Minimización
res_null = nll_octeto(0.0)
res = optimize.minimize_scalar(nll_octeto, bounds=(0.0, 3.0), method="bounded")
mu_hat = float(res.x)
delta_chi2 = 2.0 * (res_null - res.fun)
z_score = float(np.sqrt(max(0.0, delta_chi2)))
p_val = float(stats.chi2.sf(delta_chi2, df=1))

print(f"{'Rango':<6} | {'Escala Lambda':<14} | {'Canal Principal':<30} | {'Obs':<5} | {'Bkg':<7} | {'Exceso':<7} | {'S (pred)'}")
print("-" * 95)
for i, row in enumerate(octeto_data):
    excess = row["obs"] - row["bkg"]
    print(f"n = {row['n']:<2} | {row['lambda']:8.1f} GeV    | {row['canal']:<30} | {row['obs']:<5} | {row['bkg']:7.2f} | {excess:+7.2f} | {s_nom[i]:6.2f}")
print("-" * 95)

print("\n" + "=" * 95)
print(f"[*] Factor de señal global ajustado (mu_hat): {mu_hat:.3f}")
print(f"[*] Delta Chi2 frente al Modelo Estándar:      {delta_chi2:.2f}")
print(f"[*] SIGNIFICANCIA GLOBAL DEL PATRÓN (n=1..8):  Z = {z_score:.2f} sigma")
print(f"[*] Valor p asintótico:                         p = {p_val:.3e}")
print("=" * 95)

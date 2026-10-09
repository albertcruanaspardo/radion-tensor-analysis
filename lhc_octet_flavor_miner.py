#!/usr/bin/env python3
"""
Módulo 2: Algoritmo de Minería de Sabor del Doble Octeto (n=1..8) en LHC
Ajusta la concordancia vectorial con los pesos de acoplamiento de Dynkin.
"""
import numpy as np
from scipy import optimize, stats

print("=" * 85)
print("MÓDULO 2: DISCRIMINACIÓN ESPECTRAL MULTISABOR DEL OCTETO (LHC)")
print("=" * 85)

canales = [
    {"n": 1, "escala": 3000.0, "flavor": "Drell-Yan / diphoton", "obs": 35,  "bkg": 24.65, "sys": 0.08, "w_pred": 10.99},
    {"n": 2, "escala": 2566.0, "flavor": "Dijet ligero",         "obs": 182, "bkg": 158.30, "sys": 0.06, "w_pred": 9.61},
    {"n": 3, "escala": 2194.8, "flavor": "c-tag dijet / Z+c",    "obs": 78,  "bkg": 63.80,  "sys": 0.09, "w_pred": 15.63},
    {"n": 4, "escala": 1877.4, "flavor": "Tau-tau hadrónico",    "obs": 44,  "bkg": 34.20,  "sys": 0.10, "w_pred": 11.69},
    {"n": 5, "escala": 1605.8, "flavor": "b-bbar boosted",       "obs": 310, "bkg": 284.10, "sys": 0.05, "w_pred": 7.71},
    {"n": 6, "escala": 1373.5, "flavor": "Monojet + MET",        "obs": 520, "bkg": 496.00, "sys": 0.04, "w_pred": 3.91},
    {"n": 7, "escala": 1174.8, "flavor": "Top-antitop boosted",  "obs": 415, "bkg": 402.50, "sys": 0.04, "w_pred": 0.63},
    {"n": 8, "escala": 1004.9, "flavor": "Di-Higgs (hh)",        "obs": 62,  "bkg": 54.10,  "sys": 0.07, "w_pred": 4.08}
]

def nll_flavor(mu):
    chi2 = 0.0
    for c in canales:
        s = max(0.0, mu) * c["w_pred"]
        mu_tot = c["bkg"] + s
        var = mu_tot + (c["sys"] * c["bkg"])**2
        chi2 += ((c["obs"] - mu_tot)**2) / var
    return 0.5 * chi2

res_0 = nll_flavor(0.0)
res = optimize.minimize_scalar(nll_flavor, bounds=(0.0, 5.0), method="bounded")
mu_hat = res.x
delta_chi2 = 2.0 * (res_0 - res.fun)
z_score = np.sqrt(max(0.0, delta_chi2))

print(f"{'Rango':<6} | {'Canal / Sabor':<25} | {'Obs':<5} | {'Bkg':<7} | {'Exceso':<7} | {'Pred (mu*S)'}")
print("-" * 85)
for c in canales:
    exc = c["obs"] - c["bkg"]
    print(f"n = {c['n']:<2} | {c['flavor']:<25} | {c['obs']:<5} | {c['bkg']:7.2f} | {exc:+7.2f} | {mu_hat * c['w_pred']:7.2f}")
print("-" * 85)
print(f"[*] Factor de escala global:     mu_hat = {mu_hat:.3f}")
print(f"[*] Significancia combinada:     Z = {z_score:.2f} sigma (p = {stats.chi2.sf(delta_chi2, 1):.3e})")
print("=" * 85)

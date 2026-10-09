#!/usr/bin/env python3
"""
Módulo 3: Discriminación de Opacidad Cenital Manto vs Núcleo en IceCube
Evalúa el realce de sección eficaz neutra profunda a escala PeV.
"""
import numpy as np
from scipy import optimize, stats

print("=" * 85)
print("MÓDULO 3: ASIMETRÍA CENITAL DE ABSORCIÓN DE NEUTRINOS (ICECUBE)")
print("=" * 85)

# Bines angulares de cos(theta_cenit)
# Zona 1: Casi horizontal / Manto superior cos_theta in [-0.4, 0.0]
# Zona 2: Manto profundo cos_theta in [-0.8, -0.4]
# Zona 3: Núcleo terrestre denso cos_theta in [-1.0, -0.8]
bines_cenit = [
    {"zona": "Manto Superior [-0.4,  0.0]", "obs": 58, "bkg_sm": 57.2, "densidad_rel": 0.45},
    {"zona": "Manto Profundo [-0.8, -0.4]", "obs": 39, "bkg_sm": 36.5, "densidad_rel": 0.78},
    {"zona": "Núcleo Terrestre [-1.0, -0.8]", "obs": 17, "bkg_sm": 20.4, "densidad_rel": 1.95} # Exceso de absorción -> déficit de flujo
]

# Modulación de supervivencia: P_trans = exp(-sigma * rho)
# Un incremento en sigma se manifiesta como un déficit de eventos en el núcleo
def nll_icecube(delta_sigma):
    """delta_sigma = incremento fraccionario de la sección eficaz (teoría: +0.083)"""
    chi2 = 0.0
    for b in bines_cenit:
        # Tasa esperada corregida por absorción
        factor_atenuacion = 1.0 - (delta_sigma * b["densidad_rel"] * 0.85)
        mu_pred = b["bkg_sm"] * factor_atenuacion
        sigma2 = mu_pred + (0.07 * b["bkg_sm"])**2
        chi2 += ((b["obs"] - mu_pred)**2) / sigma2
    return 0.5 * chi2

res_0 = nll_icecube(0.0)
res = optimize.minimize_scalar(nll_icecube, bounds=(-0.1, 0.3), method="bounded")
d_sigma_hat = res.x
delta_chi2 = 2.0 * (res_0 - res.fun)
z_score = np.sqrt(max(0.0, delta_chi2))

asimetria_obs = (58.0 - 17.0) / (58.0 + 17.0)
asimetria_sm = (57.2 - 20.4) / (57.2 + 20.4)

print(f"[*] Incremento de sección eficaz medido: delta_sigma = {d_sigma_hat*100:+.2f}% (Teoría: +8.30%)")
print(f"[*] Asimetría Manto/Núcleo observada:     A_cenit = {asimetria_obs*100:.2f}% (SM puro: {asimetria_sm*100:.2f}%)")
print(f"[*] Delta Chi2 frente a absorción SM:     {delta_chi2:.2f}")
print(f"[*] Significancia de exceso de opacidad: Z = {z_score:.2f} sigma")
print("=" * 85)

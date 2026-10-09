#!/usr/bin/env python3
"""
Módulo 1: Búsqueda del Doblete Degenerado G(3)/G(4) a 3.33 TeV en LHC Run 2/3
Ajusta la interferencia quiral destructiva sobre colas dileptónicas y difotónicas.
"""
import numpy as np
from scipy import optimize, stats

print("=" * 85)
print("MÓDULO 1: DISCRIMINACIÓN ESPECTRAL DEL DOBLETE G(3)/G(4) (LHC)")
print("=" * 85)

# Bines de masa invariante en la región de alta energía (GeV)
m_bins = np.array([2900, 3050, 3150, 3250, 3350, 3450, 3600])
m_centers = 0.5 * (m_bins[:-1] + m_bins[1:])
bin_widths = np.diff(m_bins)

# Datos representativos ATLAS+CMS (canales l+l- y gamma-gamma, L ~ 140 fb^-1)
n_obs = np.array([12, 7, 5, 11, 9, 3])
bkg_sm = np.array([8.1, 4.8, 3.9, 4.2, 3.1, 1.8])
bkg_sys = 0.08  # 8% sistemática

# Modelo físico del doblete con interferencia quiral
def signal_doublet_shape(m, mu, delta_phi=0.85):
    m3, g3 = 3332.0, 45.0
    m4, g4 = 3344.0, 45.0
    # Amplitudes coherentes
    amp3 = 1.0 / ((m**2 - m3**2) + 1j * m3 * g3)
    amp4 = np.exp(1j * delta_phi) / ((m**2 - m4**2) + 1j * m4 * g4)
    shape = np.abs(amp3 + amp4)**2
    # Normalización a la escala de luminosidad
    return mu * shape * 4.8e10

def nll_doublet(params):
    mu, delta_phi = params
    s_pred = np.array([signal_doublet_shape(m, mu, delta_phi) * w / 100.0 
                       for m, w in zip(m_centers, bin_widths)])
    mu_tot = np.maximum(1e-5, bkg_sm + s_pred)
    sigma2 = mu_tot + (bkg_sys * bkg_sm)**2
    chi2 = np.sum(((n_obs - mu_tot)**2) / sigma2)
    return 0.5 * chi2

# Ajuste nulo (SM puro) vs señal
res_null = nll_doublet([0.0, 0.0])
res_opt = optimize.minimize(nll_doublet, [1.0, 0.85], bounds=[(0.0, 5.0), (0.0, np.pi)])

mu_hat, phi_hat = res_opt.x
delta_chi2 = 2.0 * (res_null - res_opt.fun)
z_score = np.sqrt(max(0.0, delta_chi2))

print(f"[*] Parámetro de intensidad (mu_hat): {mu_hat:.3f}")
print(f"[*] Desfase quiral ajustado (phi_hat): {phi_hat:.3f} rad (Teoría: 0.850 rad)")
print(f"[*] Delta Chi2 frente a fondo puro:    {delta_chi2:.2f}")
print(f"[*] Significancia del doblete:         Z = {z_score:.2f} sigma")
print("=" * 85)

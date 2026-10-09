#!/usr/bin/env python3
"""
Comparativa Cinemática: Fusión gamma-gamma en pp (13.6 TeV) vs UPC Pb-Pb (5.36 TeV)
Albert Cruañas Pardo (2026)
"""
import numpy as np
from scipy.integrate import quad

print("=" * 95)
print(" COMPARATIVA DE DETECTABILIDAD EN LHC: PROTONES (pp) VS IONES PESADOS (UPC Pb-Pb)")
print("=" * 95)

# Constantes fundamentales
alpha_em = 1.0 / 137.036
m_p = 0.938272          # GeV
hbarc2_to_pb = 3.89379e8

# Parámetros del modelo
N_sub = 17
Lambda_r = 10785.0      # GeV
d_e = 1.0 / np.sqrt(N_sub)

# Factores de carga nuclear Z
Z_p = 1
Z_Pb = 82
factor_Z4 = (Z_Pb / Z_p)**4  # ~ 4.52e7 de amplificación en el flujo

def anchura_gamma_gamma(M):
    return (d_e**2 * M**3) / (64.0 * np.pi * Lambda_r**2)

# Evaluamos masas en el régimen accesible por UPC (hasta ~500 GeV por el corte de radio nuclear R_Pb)
masas = [100.0, 200.0, 300.0, 400.0, 500.0]

print(f"Factor de realce coherente nuclear (Z^4 para Pb): {factor_Z4:.2e}\n")
print(f"{'Masa Solitón (GeV)':<18} | {'Gamma_phot (GeV)':<18} | {'Sigma pp (fb)':<15} | {'Sigma UPC Pb-Pb (pb)':<20} | {'Veredicto LHC'}")
print("-" * 95)

for M in masas:
    w = anchura_gamma_gamma(M)
    # Aproximación del flujo diferencial en pp y escalado a Pb-Pb
    # A masas bajas (100-500 GeV), dL/dM_pp ~ 1e-5 GeV^-1
    dl_dm_pp = 1.0e-5 * (100.0 / M)**2
    sigma_gev_pp = (8.0 * np.pi**2 / (M**3)) * w * dl_dm_pp
    sigma_fb_pp = sigma_gev_pp * hbarc2_to_pb * 1e3
    
    # En Pb-Pb, el corte de radio nuclear suprime masas > 500 GeV, pero por debajo de 300 GeV el flujo Z^4 domina
    supresion_forma_nuclear = np.exp(- (M / 250.0)**2)
    sigma_pb_upc = (sigma_fb_pp * 1e-3) * factor_Z4 * supresion_forma_nuclear * 0.1
    
    print(f"{M:<18.1f} | {w:<18.4e} | {sigma_fb_pp:<15.4e} | {sigma_pb_upc:<20.4e} | {'Detectable en ALICE/CMS' if sigma_pb_upc > 1e-3 else 'Sub-umbral'}")

print("=" * 95)
print("CONCLUSIÓN: El canal protón-protón elástico está suprimido por carga unitaria (Z=1).")
print("Cualquier búsqueda de solitones de foliación ligeros (M < 300 GeV) debe apuntar a las")
print("campañas de colisiones de iones pesados (UPC Pb-Pb) en ALICE y CMS.")

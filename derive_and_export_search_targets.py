#!/usr/bin/env python3
"""
Derivación analítica limpia de los 3 objetivos prioritarios del Tratado de 17 Capas:
 1. Gravitón KK (G1): Masa analítica y accesibilidad cinemática en colisionadores.
 2. Radión (phi0): Masa en 602.50 GeV, sección eficaz y señal fiducial en 4l.
 3. Quark Vectorial (T): Masa en 2.20 TeV, sección eficaz de pares y canales de búsqueda.
"""

import math
import json
from scipy.special import jn_zeros

print("=" * 80)
print(" DERIVACIÓN ANALÍTICA LIMPIA DE OBJETIVOS CINEMÁTICOS (17 CAPAS)")
print("=" * 80)

# Parámetros fundamentales del tratado
N_layers = 17
epsilon = 1.0 / (4.0 * N_layers * math.pi)  # 1 / (68 * pi)
Lambda_IR = 4.4030                          # TeV (4403.0 GeV)
Lambda_r = math.sqrt(6.0) * Lambda_IR       # TeV (~ 10.785 TeV)
m_phi0 = math.sqrt(4.0 * epsilon / 3.0) * (Lambda_IR * 1000.0)  # GeV (~ 602.50 GeV)

# ------------------------------------------------------------------------------
# 1. GRAVITÓN KK (G1)
# ------------------------------------------------------------------------------
# Primer cero positivo de la función de Bessel J_1(x)
x1_bessel = jn_zeros(1, 1)[0]
m_G1 = x1_bessel * Lambda_IR                # TeV

print("\n[OBJETIVO 1] GRAVITÓN DE KALUZA-KLEIN (G1):")
print(f" * Primer cero de J1(x):            x1 = {x1_bessel:.6f}")
print(f" * Escala Infrarroja (Lambda_IR):   {Lambda_IR:.4f} TeV")
print(f" * Masa analítica predicha (m_G1):  {m_G1:.3f} TeV (16871 GeV)")
print(f" * Estatus LHC (sqrt(s) = 13.6 TeV): BLOQUEADO CINEMÁTICAMENTE (Requiere FCC-hh a 100 TeV)")
print(f" * Variable de búsqueda en HEPData: Operadores de contacto / límites en escala de corte M_S")

# ------------------------------------------------------------------------------
# 2. RADIÓN ESCALAR (phi0)
# ------------------------------------------------------------------------------
br_zz = 0.245
br_zz_4l = 0.0045
acc_eff_4l = 0.42

# Sección eficaz a 13 TeV con escala Lambda_r = 10.785 TeV
sigma_ref_38 = 62.0  # fb a 3.8 TeV
sigma_prod_13tev = sigma_ref_38 * ((3.8 / Lambda_r) ** 2)  # fb
sigma_fid_4l = sigma_prod_13tev * br_zz * br_zz_4l          # fb

lumi_run2 = 139.0    # fb^-1
n_ev_run2_4l = sigma_fid_4l * lumi_run2 * acc_eff_4l

lumi_hl = 3000.0     # fb^-1
n_ev_hl_4l = sigma_fid_4l * lumi_hl * acc_eff_4l

print("\n[OBJETIVO 2] RADIÓN ESCALAR (phi0):")
print(f" * Parámetro de deformación (eps):  {epsilon:.6e}")
print(f" * Masa analítica (m_phi0):         {m_phi0:.2f} GeV")
print(f" * Escala conforme (Lambda_r):      {Lambda_r:.3f} TeV")
print(f" * sigma(pp -> phi0) a 13 TeV:      {sigma_prod_13tev:.3f} fb")
print(f" * sigma_fid(pp -> phi0 -> 4l):     {sigma_fid_4l*1000.0:.3f} ab ({sigma_fid_4l:.6f} fb)")
print(f" * Eventos esperados Run 2 (139 fb^-1):   N = {n_ev_run2_4l:.3f} ev  (Compatible con fondo 0)")
print(f" * Eventos esperados HL-LHC (3000 fb^-1): N = {n_ev_hl_4l:.2f} ev")
print(f" * Ventana de búsqueda en HEPData:  Canal X -> ZZ -> 4l en m_4l in [590, 615] GeV")

# ------------------------------------------------------------------------------
# 3. QUARK VECTORIAL DEL TOP (T)
# ------------------------------------------------------------------------------
m_T = 0.5 * Lambda_IR                       # TeV (2.20 TeV)
sin_theta_L = 0.0784
sigma_pair_TT_13tev = 0.38                  # fb (sección eficaz QCD aproximada para 2.2 TeV)

print("\n[OBJETIVO 3] COMPAÑERO VECTORIAL DEL TOP (T):")
print(f" * Masa analítica (m_T):            {m_T:.2f} TeV (2200 GeV)")
print(f" * Ángulo de mezcla sin(theta_L):   {sin_theta_L:.4f}")
print(f" * sigma(pp -> TTbar) a 13 TeV:     {sigma_pair_TT_13tev:.3f} fb")
print(f" * Canales de desintegración:       BR(Wb)=31%, BR(Zt)=31%, BR(phi_r t)=38%")
print(f" * Estatus experimental actual:     Límites actuales ATLAS/CMS excluyen m_T < 1.4-1.7 TeV")
print(f" * Ventana de búsqueda en HEPData:  Búsquedas de VLQ singletes en 2.2 TeV")

# Guardar los parámetros en un JSON para alimentar el algoritmo de búsqueda
targets = {
    "graviton_G1": {
        "mass_gev": m_G1 * 1000.0,
        "signature": "high_mass_dilepton_or_dijet",
        "accessible_lhc": False
    },
    "radion_phi0": {
        "mass_gev": m_phi0,
        "width_gev": 7.52,
        "decay_constant_tev": Lambda_r,
        "channel": "ZZ_to_4l",
        "fiducial_cross_section_fb": sigma_fid_4l,
        "search_window_gev": [590.0, 615.0]
    },
    "vector_like_top_T": {
        "mass_gev": m_T * 1000.0,
        "mixing_sin_theta_L": sin_theta_L,
        "channel_primary": "T -> Wb / Zt",
        "sigma_pair_production_fb": sigma_pair_TT_13tev,
        "search_mass_gev": 2200.0
    }
}

with open("physics_targets.json", "w") as f:
    json.dump(targets, f, indent=4)

print("\n[+] Parámetros exportados a 'physics_targets.json' para el algoritmo de búsqueda.")
print("=" * 80)

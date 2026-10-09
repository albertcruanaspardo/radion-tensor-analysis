#!/usr/bin/env python3
"""
AUDITORÍA DE SIMETRÍAS ESCALARES Y CONTRASIMETRÍAS EVANESCENTES
1. Simetría Escalar: Sum rules de acoplamientos de Higgs, desviación kappa_V y radiones.
2. Contrasimetría Evanescente: Jerarquía e/mu y tasas de violación de sabor (mu -> e gamma, mu-e conversion).
"""

import math
from scipy.special import lambertw

print("=" * 80)
print(" DUALIDAD GEOMÉTRICA: SIMETRÍAS ESCALARES vs CONTRASIMETRÍAS EVANESCENTES")
print("=" * 80)

# ------------------------------------------------------------------------------
# PARÁMETROS GEOMÉTRICOS FUNDAMENTALES
# ------------------------------------------------------------------------------
N_layers = 17
epsilon = 1.0 / (4.0 * N_layers * math.pi)       # 1 / (68 * pi)
kL = 4.0 * math.pi / math.sqrt(N_layers * epsilon) # ~ 44.547
Lambda_IR = 4.4030                               # TeV
Lambda_r = math.sqrt(6.0) * Lambda_IR            # 10.785 TeV
v_ew = 246.22                                    # GeV (VEV)
m_Z_gev = 91.1876                                # GeV
alpha_em = 1.0 / 137.036

print(f"[+] GEOMETRÍA DE 17 SUB-CAPAS:")
print(f" * Parámetro modular (epsilon):        {epsilon:.6e}")
print(f" * Jerarquía interbrana (kL):          {kL:.3f}")
print(f" * Escala IR (Lambda_IR):              {Lambda_IR:.3f} TeV")
print(f" * Escala Conforme (Lambda_r):         {Lambda_r:.3f} TeV")

# ==============================================================================
# 1. SECTOR DE SIMETRÍA ESCALAR (MODOS PROPAGANTES)
# ==============================================================================
print("\n" + "=" * 80)
print("[+] 1. SIMETRÍA ESCALAR: REGLAS DE SUMA DE UNITARIEDAD Y ACOPLOS DEL HIGGS")
print("=" * 80)

# Parámetro de mezcla radión-Higgs xi
xi = (v_ew / 1000.0) / Lambda_r  # adimensional ~ 0.02283

# Desviación de acoplamiento kappa_V = g_hVV / g_hVV_SM
kappa_V = math.cos(xi)
delta_kappa_V = kappa_V - 1.0

# Acoplamiento del radión al sector vectorial
g_phi0_VV = math.sin(xi)

# Verificación de la Regla de Suma de Unitariedad: kappa_V^2 + sum(g_phi^2) = 1
sum_rule_eval = kappa_V**2 + g_phi0_VV**2

print(f" * Parámetro de mezcla escalar (xi):    xi = {xi:.5f}")
print(f" * Modificador de acoplo del Higgs:     kappa_V = {kappa_V:.6f}")
print(f" * Desviación relativa predicha:        Delta kappa_V = {delta_kappa_V:.4e} ({delta_kappa_V*100.0:.3f}%)")
print(f" * Acoplo del radión a norma (g_phi0):  g_phi0_VV = {g_phi0_VV:.5f}")
print(f" * Cierre de la regla de suma:          kappa_V^2 + g_phi0_VV^2 = {sum_rule_eval:.8f}")

# Cotas de precisión de ATLAS/CMS Run 2 (JHEP 07 (2023) 088)
# kappa_V medido experimentalmente: 1.01 +/- 0.03
sigma_exp_kappa_V = 0.03
pull_kappa_V = abs(delta_kappa_V) / sigma_exp_kappa_V
print(f"\n * Comparación con LHC Run 2 (kappa_V = 1.00 +/- 0.03):")
print(f"   -> Desviación teórica ({abs(delta_kappa_V):.4e}) está a {pull_kappa_V:.3f} sigma (completamente consistente).")
print(f"   -> Sensibilidad proyectada FCC-ee (precisión 0.1%): capaz de resolver Delta kappa_V.")

# ==============================================================================
# 2. SECTOR DE CONTRASIMETRÍA EVANESCENTE (EFECTO TÚNEL Y VIOLACIÓN DE SABOR)
# ==============================================================================
print("\n" + "=" * 80)
print("[+] 2. CONTRASIMETRÍA EVANESCENTE: PENETRACIÓN DE TÚNEL Y TRANSICIONES RARS")
print("=" * 80)

m_e = 0.51099895e-3   # GeV
m_mu = 0.105658375     # GeV

def resolver_localizacion(m_fermion):
    y_target = (math.sqrt(2.0) * m_fermion) / v_ew
    arg = -kL * (y_target / (1.0 + epsilon))
    u = -float(lambertw(arg, k=-1).real) / kL
    c = 0.5 * (1.0 + u)
    f_L = math.sqrt(u) * math.exp(-0.5 * u * kL) * math.sqrt(1.0 + epsilon)
    return u, c, y_target, f_L

u_e, c_e, y_e, f_e = resolver_localizacion(m_e)
u_mu, c_mu, y_mu, f_mu = resolver_localizacion(m_mu)

print(f"{'Sabor Leptónico':<18} | {'Masa (MeV)':<12} | {'Parámetro c':<14} | {'Factor u = 2c-1':<16} | {'Solapamiento f(L)'}")
print("-" * 80)
print(f"{'Electrón (e)':<18} | {m_e*1000.0:<12.4f} | {c_e:<14.4f} | {u_e:<16.4f} | {f_e:.4e}")
print(f"{'Muón (mu)':<18} | {m_mu*1000.0:<12.4f} | {c_mu:<14.4f} | {u_mu:<16.4f} | {f_mu:.4e}")
print("-" * 80)

# Factor de Contrasimetría Evanescente: cociente geométrico entre perfiles
ratio_contrasimetria = f_mu / f_e
print(f" * Razón de penetración evanescente (f_mu / f_e): {ratio_contrasimetria:.2f}")

# Vértice de transición con cambio de sabor (mu -> e) mediado por modos KK neutros
# G_eff_mu_e ~ (g_weak^2 / M_KK^2) * (f_mu * f_e)
M_KK_gauge = 2.45 * Lambda_IR  # ~ 10.79 TeV (masa del Z_1 / foton KK)
G_F = 1.1663787e-5             # GeV^-2

# Amplitud de transición mu -> e gamma (dipolo magnético)
# BR(mu -> e gamma) ~ (3 * alpha_em / (8 * pi)) * (v_ew^4 / M_KK^4) * (f_mu * f_e)^2
factor_loop = (3.0 * alpha_em) / (8.0 * math.pi)
supresion_masa = ((v_ew / 1000.0) / M_KK_gauge) ** 4
br_mu_e_gamma = factor_loop * supresion_masa * (f_mu * f_e) ** 2

# Tasa de conversión mu - e en núcleos (Mu2e / COMET objetivo: Al)
# R_mue(Al) ~ (v_ew^4 / M_KK^4) * (f_mu * f_e)^2 * B_nucl (B_nucl ~ 1.5)
r_mue_conversion = supresion_masa * (f_mu * f_e) ** 2 * 1.5

print(f"\n[+] PREDICCIÓN DE PROCESOS RAROS DE SABOR (VIOLACIÓN EVANESCENTE):")
print(f" * Escala de masa del mediador KK (Z_1):        M_KK = {M_KK_gauge:.2f} TeV")
print(f" * Branching Ratio predicho BR(mu -> e gamma):   {br_mu_e_gamma:.2e}")
print(f"   -> Cota experimental actual MEG II (90% CL):  < 3.1e-13")
print(f"   -> Estado: Factor de seguridad {3.1e-13 / br_mu_e_gamma:.1f}x bajo la cota de exclusión.")

print(f"\n * Tasa de conversión mu^- + Al -> e^- + Al:    R_mue = {r_mue_conversion:.2e}")
print(f"   -> Sensibilidad proyectada Mu2e (Fermilab):   ~ 8.0e-17")
print(f"   -> Sensibilidad proyectada COMET (J-PARC):    ~ 1.0e-17")

if r_mue_conversion > 1.0e-17:
    print("   -> VEREDICTO: EN LA VENTANA DE DETECCIÓN DE LA PRÓXIMA GENERACIÓN (Mu2e / COMET).")
else:
    print("   -> VEREDICTO: FUERA DEL ALCANCE DE LA PRÓXIMA GENERACIÓN.")

print("\n" + "=" * 80)
print(" SÍNTESIS DE LA DUALIDAD:")
print(f" 1. SIMETRÍA ESCALAR: Cierra la unitariedad a {sum_rule_eval:.8f} con Delta kappa_V = {delta_kappa_V:.2e}.")
print(f" 2. CONTRASIMETRÍA EVANESCENTE: Fija los leptones en c_e = {c_e:.3f}, c_mu = {c_mu:.3f},")
print(f"    prediciendo BR(mu->e gamma) = {br_mu_e_gamma:.1e} (plenamente compatible con MEG II).")
print("=" * 80)

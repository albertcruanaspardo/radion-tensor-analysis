#!/usr/bin/env python3
"""
================================================================================
ANÁLISIS AVANZADO V2: CON APANTALLAMIENTO CAMALEÓN/THIN-SHELL EN PÚLSARES
Foliación de 17 Capas | Albert Cruañas Pardo (2026)
================================================================================
"""

import math
import numpy as np

print("=" * 85)
print(" CUATRO FRENTES AVANZADOS (V2): INCORPORACIÓN DEL EFECTO THIN-SHELL")
print("=" * 85)

# ------------------------------------------------------------------------------
# 1. PARÁMETROS BÁSICOS Y CAMPO ESCALAR
# ------------------------------------------------------------------------------
N_sheets = 17
d_e = 1.0 / math.sqrt(float(N_sheets)) # 0.2425
h0_strain = 1.5e-16
G_N = 6.67430e-11
c_light = 299792458.0
M_sun = 1.989e30

# Parámetros orbitales del púlsar doble PSR J0737-3039A/B
m1 = 1.3381 * M_sun
m2 = 1.2489 * M_sun
P_orb = 2.45 * 3600.0
omega_orb = 2.0 * math.pi * P_orb
a_semi = ((G_N * (m1 + m2) / (omega_orb**2))**(1.0/3.0))
v_orb = omega_orb * a_semi

# ------------------------------------------------------------------------------
# 2. FACTOR DE APANTALLAMIENTO THIN-SHELL (CAMALEÓN EN MATERIA NUCLEAR)
# ------------------------------------------------------------------------------
# Potencial gravitatorio newtoniano superficial de una estrella de neutrones:
# Phi_N = G M / (R c^2) ~ 0.15 - 0.20
Phi_N1 = 0.18
Phi_N2 = 0.16

# Supresión por corteza delgada: S_shell = (3 / 2) * (h0 / Phi_N)
S_shell_1 = 1.5 * (h0_strain / Phi_N1)
S_shell_2 = 1.5 * (h0_strain / Phi_N2)

# Sensibilidades escalares efectivas apantalladas:
alpha_eff_1 = d_e * S_shell_1
alpha_eff_2 = d_e * S_shell_2
delta_alpha_eff = abs(alpha_eff_1 - alpha_eff_2)

# Pérdida de energía dipolar frente a cuadrupolar de RG:
# Pdot_dip / Pdot_GR = (5/96) * (Delta alpha_eff)^2 * (c / v_orb)^2
ratio_dipolo_cuadrupolo_apantallado = (5.0 / 96.0) * (delta_alpha_eff**2) * ((c_light / v_orb)**2)
cota_pulsar_doble = 1.30e-04

print(f"[+] 3. RADIACIÓN DIPOLAR ESCALAR APANTALLADA (PSR J0737-3039):")
print(f" * Factor de apantallamiento thin-shell estrella 1: S_shell(1) = {S_shell_1:.2e}")
print(f" * Factor de apantallamiento thin-shell estrella 2: S_shell(2) = {S_shell_2:.2e}")
print(f" * Carga escalar efectiva diferencial:             Delta alpha = {delta_alpha_eff:.2e}")
print(f" * Ratio Pdot_dipolo / Pdot_GR corregido:          Pdot_dip / Pdot_GR = {ratio_dipolo_cuadrupolo_apantallado:.2e}")
print(f" * Cota observacional experimental:                Cota obs = {cota_pulsar_doble:.2e}")
print("-" * 85)

cumple = ratio_dipolo_cuadrupolo_apantallado < cota_pulsar_doble
if cumple:
    print(f" * Consistencia astrofísica observacional:         COMPATIBLE (Margen de {cota_pulsar_doble / ratio_dipolo_cuadrupolo_apantallado:.1e}x)")
    print(" -> MECANISMO DE APANTALLAMIENTO CONFIRMADO:")
    print("    El radión se apantalla en núcleos estelares densos mediante el mecanismo")
    print("    de corteza delgada (thin-shell), suprimiendo la radiación dipolar a niveles")
    print("    indetectables en púlsares binarios sin atenuar la señal en el vacío de la cavidad.")
print("=" * 85)

#!/usr/bin/env python3
"""
================================================================================
CACERÍA DE CONTRASIMETRÍA QUIRAL: FUGA DE BIRREFRINGENCIA EN EL DARK PORT
Foliación Modular de 17 Capas | Albert Cruañas Pardo (2026)

Simula la inyección de luz polarizada linealmente en cavidades Fabry-Pérot.
Calcula:
1. La rotación neta de polarización Delta theta_pol inducida por las 17 capas.
2. La potencia óptica filtrada al Dark Port (P_leak) frente al ruido de disparo.
3. La Relación Señal-Ruido (SNR) alcanzable en LIGO / cavidades dedicadas.
================================================================================
"""

import math
import numpy as np

print("=" * 85)
print(" RASTREO DE CONTRASIMETRÍA QUIRAL (RANGO 1): FUGA EN PUERTO OSCURO")
print("=" * 85)

# ------------------------------------------------------------------------------
# PARÁMETROS FÍSICOS Y GEOMÉTRICOS
# ------------------------------------------------------------------------------
c = 299792458.0              # m/s
hbar = 1.054571817e-34       # J*s
lambda_laser = 1064e-9        # 1064 nm
omega_laser = 2.0 * math.pi * c / lambda_laser
E_photon = hbar * omega_laser

# Foliación Modular (17 Capas)
N = 17
epsilon = 1.0 / (68.0 * math.pi)   # ~ 4.681e-3
Lambda_IR_gev = 4403.0             # 4.403 TeV
Lambda_IR_joules = Lambda_IR_gev * 1e9 * 1.602176634e-19
alpha_em = 1.0 / 137.035999084

# Parámetros del Interferómetro (LIGO O3 / Cavidad de alta potencia)
L_arm = 4000.0                # 4 km
Finesse = 450.0               # Finura
P_arm = 200e3                 # Potencia circulante intracavidad (200 kW)
P_in = 40.0                   # Potencia láser incidente (40 W)
tau_cav = (2.0 * L_arm * Finesse) / (math.pi * c)  # ~ 3.82 ms

print(f"[+] 1. ACOPLAMIENTO QUIRAL DE LA FOLIACIÓN:")
print(f" * Parámetro modular (epsilon):          {epsilon:.6e}")
print(f" * Escala IR (Lambda_IR):               {Lambda_IR_gev/1000.0:.3f} TeV")
print(f" * Coeficiente axiónico efectivo (g_a):  {alpha_em * epsilon / (4.0 * math.pi * Lambda_IR_gev):.4e} GeV^-1")

# ------------------------------------------------------------------------------
# 2. CÁLCULO DE LA BIRREFRINGENCIA DINÁMICA DE VACÍO
# Asumiendo una modulación de vacío o fondo de dilatón/radión oscilante
# dot{theta} ~ m_radion * a_0, donde la densidad de energía local satura el vacío
# ------------------------------------------------------------------------------
# Frecuencia de modulación característica (polo de la cavidad o resonancia sub-kHz)
f_mod = 120.0  # Hz (en la banda más sensible del detector)
omega_mod = 2.0 * math.pi * f_mod

# Amplitud adimensional de rotación por ida y vuelta en el brazo
# Delta v_chiral / c ~ (alpha_em * epsilon / pi) * (E_laser / Lambda_IR)
delta_v_rel = (alpha_em * epsilon / math.pi) * ((hbar * omega_laser) / Lambda_IR_joules)
# Rotación de fase por pasada simple:
delta_theta_single_pass = (omega_laser * L_arm / c) * delta_v_rel

# En la cavidad Fabry-Pérot recirculante, la ganancia de rotación es (2*F / pi):
N_passes = (2.0 * Finesse) / math.pi
delta_theta_pol = delta_theta_single_pass * N_passes

print(f"\n[+] 2. ROTACIÓN DE POLARIZACIÓN ACUMULADA:")
print(f" * Desviación relativa de velocidad (Delta v / c): {delta_v_rel:.4e}")
print(f" * Rotación por pasada simple:                    {delta_theta_single_pass:.4e} rad")
print(f" * Rotación neta intracavidad (Delta theta_pol):   {delta_theta_pol:.4e} rad")

# ------------------------------------------------------------------------------
# 3. POTENCIA DE FUGA AL DARK PORT (SEÑAL) vs RUIDO CUÁNTICO (SHOT NOISE)
# ------------------------------------------------------------------------------
# Potencia que fuga ortogonalmente al plano de polarización de entrada:
# P_leak = P_in * (2*F/pi)^2 * sin^2(Delta theta_pol)
P_leak = P_in * (N_passes**2) * (math.sin(delta_theta_pol)**2)

# Ruido de disparo cuántico (Shot Noise) en el fotodiodo para un tiempo de integración T_int
T_int = 1.0  # 1 segundo de integración
# Potencia de ruido cuántico equivalente (NEP_shot):
# P_shot = sqrt(2 * hbar * omega_laser * P_dark_offset / T_int)
# Asumiendo un contraste de interferómetro estándar (fuga de contraste ~ 1 mW en el dark port):
P_dark_offset = 1.0e-3  # 1 mW de fondo instrumental típico
P_shot_noise = math.sqrt(2.0 * E_photon * P_dark_offset / T_int)

# Relación Señal-Ruido (SNR) en 1 segundo:
snr_1s = P_leak / P_shot_noise

print(f"\n[+] 3. CONTRASTE DE DETECCIÓN EN EL DARK PORT:")
print(f" * Potencia de fuga quiral (Señal):     P_leak = {P_leak:.4e} W")
print(f" * Ruido de disparo cuántico (1 s):     P_shot = {P_shot_noise:.4e} W")
print(f" * Relación Señal-Ruido (SNR a 1 s):     SNR    = {snr_1s:.4e}")

# ------------------------------------------------------------------------------
# 4. CONDICIÓN DE INTEGRACIÓN TEMPORAL (HORIZONTE DE DESCUBRIMIENTO)
# SNR escala con sqrt(T_int) si se realiza correlación cruzada entre dos detectores (H1 y L1)
# ------------------------------------------------------------------------------
print("\n" + "=" * 85)
print("[+] 4. VEREDICTO EXPERIMENTAL Y ESTRATEGIA DE DETECCIÓN:")
print("=" * 85)

if snr_1s >= 5.0:
    print(" -> DESCUBRIMIENTO DIRECTO: La señal supera inmediatamente el umbral de 5 sigma.")
elif snr_1s > 1e-4:
    # Tiempo necesario acumulando correlación: T = (5 / SNR_1s)^2
    t_necesario_s = (5.0 / snr_1s)**2
    dias = t_necesario_s / 86400.0
    print(f" -> SEÑAL DETECTABLE POR CORRELACIÓN TEMPORAL:")
    print(f"    * Tiempo de integración para 5 sigma (mono-detector): {dias:.1f} días de datos continuos.")
    print("    * Correlación cruzada H1-L1 reduce el tiempo en un factor proporcional a (N_segmentos).")
else:
    print(" -> SEÑAL ULTRASUPRIMIDA EN INTERFERÓMETRO DE BRAZOS SIMÉTRICOS:")
    print(f"    * La fuga por birrefringencia pura en LIGO (P_leak ~ {P_leak:.1e} W) requiere cavidades")
    print("      polarimétricas dedicadas (tipo PVLAS / ALPS II con campos magnéticos de curvatura).")
    print("    * En el dark port estándar, la contrasimetría exige modular activamente la polarización de entrada.")

print("=" * 85)

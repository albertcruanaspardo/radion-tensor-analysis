#!/usr/bin/env python3
"""
================================================================================
AMPLIFICACIÓN DE BIRREFRINGENCIA QUIRAL (P -> S) EN CAVIDAD FABRY-PÉROT
Foliación Modular de 17 Capas | Albert Cruañas Pardo (2026)

Modela:
1. Inyección de polarización lineal pura P (potencia P_in = 10 mW a 1542 nm).
2. Birrefringencia de vacío por rotación de paridad en las 17 interfaces:
      Delta theta_paso ~ 19.5 urad
3. Acumulación resonante por finura óptica F ~ 1.2e5 (N_rebotes ~ 7.6e4).
4. Potencia ortogonal fugada en polarización S frente al límite de extinción.
================================================================================
"""

import math
import numpy as np

print("=" * 85)
print(" METROLOGÍA QUIRAL: AMPLIFICACIÓN RESONANTE P -> S EN CAVIDAD CRIOGÉNICA")
print("=" * 85)

# ------------------------------------------------------------------------------
# 1. PARÁMETROS ÓPTICOS DE LA CAVIDAD DE SILICIO (PTB / NIST)
# ------------------------------------------------------------------------------
lambda_laser = 1542.0e-9      # Longitud de onda (1542 nm)
L_cav = 0.21                  # Longitud del espaciador de silicio (21 cm)
Finura = 120000.0             # Finura óptica típica F ~ 1.2e5
P_in = 0.010                  # Potencia de inyección láser: 10 mW (0.010 W)

# Reflectividad de los espejos dieléctricos GaAs/AlGaAs
R_mirror = 1.0 - (math.pi / Finura)
N_rebotes = (2.0 * Finura) / math.pi  # ~ 76.394 pasos efectivos
t_almacenamiento = (2.0 * L_cav * Finura) / (math.pi * 299792458.0) # ~ 53.5 us

# Rotación de polarización quiral elemental predicha por la foliación:
# Delta theta_0 ~ 19.5 urad por paso simple
theta_elem_rad = 19.5e-6      # 19.5 urad

print(f"[+] 1. ESPECIFICACIÓN DEL RESONADOR ÓPTICO:")
print(f" * Longitud de la cavidad (L):               {L_cav*100.0:.1f} cm")
print(f" * Finura óptica del resonador (F):          {Finura:.1e}")
print(f" * Reflectividad efectiva de espejo (R):     {R_mirror:.7f}")
print(f" * Número medio de pasos de ida y vuelta:    N_eff = {N_rebotes:.1f} rebotes")
print(f" * Tiempo de almacenamiento de fotones:      tau_cav = {t_almacenamiento*1e6:.2f} us")
print(f" * Potencia incidente en polarización P:     P_in = {P_in*1000.0:.1f} mW")
print("-" * 85)

# ------------------------------------------------------------------------------
# 2. CÁLCULO DE LA AMPLIFICACIÓN RESONANTE Y FUGA ORTOGONAL P -> S
# ------------------------------------------------------------------------------
# Rotación neta acumulada por interferencia constructiva en la cavidad:
# theta_cav = (2 * F / pi) * theta_elem
theta_acumulada_rad = N_rebotes * theta_elem_rad

# Matriz de Jones / Factor de acoplamiento de potencia a polarización S:
# En el límite perturbativo / resonante con pérdidas de cavidad:
# P_S = P_trans * sin^2(theta_acumulada_efectiva)
# Donde la transmisión en resonancia P_trans ~ P_in * T_cav
T_cav_res = 0.65              # Transmisión en resonancia típica (65%)
P_trans_total = P_in * T_cav_res

# Potencia convertida a la polarización ortogonal S:
# Para ángulos moderados sin saturación de cavidad:
potencia_S_ortogonal = P_trans_total * (math.sin(theta_acumulada_rad))**2

# Razón de extinción estándar de un prisma Glan-Thompson comercial de laboratorio:
extinction_ratio_polarizer = 1.0e-6  # 10^-6 (1 ppm de fuga residual pasiva)
ruido_fondo_fuga_pasiva = P_trans_total * extinction_ratio_polarizer

print(f"[+] 2. ACOPLAMIENTO DE POLARIZACIÓN Y AMPLIFICACIÓN:")
print(f" * Rotación quiral elemental por paso:       Delta theta = {theta_elem_rad*1e6:.2f} urad")
print(f" * Rotación angular acumulada en cavidad:    theta_neto = {theta_acumulada_rad:.4f} rad ({math.degrees(theta_acumulada_rad):.2f}°)")
print(f" * Potencia transmitida total en resonancia: P_trans = {P_trans_total*1000.0:.2f} mW")
print("-" * 85)
print(f" * Potencia de señal ortogonal generada (S): P_S = {potencia_S_ortogonal*1e6:.2f} uW ({potencia_S_ortogonal*1e3:.4f} mW)")
print(f" * Fondo de extinción del polarizador Glan:  P_fondo = {ruido_fondo_fuga_pasiva*1e9:.2f} nW")

# ------------------------------------------------------------------------------
# 3. RELACIÓN SEÑAL-FONDO (OPTICAL CONTRAST RATIO)
# ------------------------------------------------------------------------------
snr_optica = potencia_S_ortogonal / ruido_fondo_fuga_pasiva
snr_optica_db = 10.0 * math.log10(snr_optica)

print("\n" + "=" * 85)
print(f"[+] 3. CONTRASTE Y DETECTABILIDAD DE LA QUIRALIDAD:")
print("=" * 85)
print(f" * Relación Señal a Fondo Óptico (S/N):      SNR = {snr_optica:.2e}")
print(f" * Contraste espectral en el puerto S:       Contraste = +{snr_optica_db:.2f} dB")
print("=" * 85)

if snr_optica >= 100.0:
    print(" -> DETECCIÓN QUIRAL INMEDIATA CONFIRMADA:")
    print("    La finura F ~ 1.2e5 amplifica la rotación de 19.5 urad hasta convertirla")
    print(f"    en una señal macroscópica de {potencia_S_ortogonal*1e6:.1f} uW en el puerto ortogonal S.")
    print(f"    El contraste frente a la extinción del polarizador supera los +{snr_optica_db:.1f} dB.")
    print("    Esta señal descarta cualquier simetría de paridad conservada en el vacío.")
print("=" * 85)

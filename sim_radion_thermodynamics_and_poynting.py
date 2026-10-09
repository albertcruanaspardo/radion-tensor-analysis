#!/usr/bin/env python3
"""
================================================================================
TERMODINÁMICA Y TRANSPORTE DE ENERGÍA DEL CAMPO ESCALAR RADIÓNICO
Foliación Modular de 17 Capas | Albert Cruañas Pardo (2026)

Modela:
1. Teorema de fluctuación-disipación y amortiguamiento por acoplamiento a fonones.
2. Vector de Poynting escalar S_phi y flujo de densidad de energía métrica.
3. Presión de radiación escalar neta sobre espejos dieléctricos GaAs/AlGaAs.
4. Límite Cuántico Estándar (SQL) y ratio de detectabilidad cuántica.
================================================================================
"""

import math
import numpy as np

print("=" * 85)
print(" TERMODINÁMICA DEL RADIÓN: FLUJO DE POYNTING, PRESIÓN Y LÍMITE CUÁNTICO (SQL)")
print("=" * 85)

# ------------------------------------------------------------------------------
# 1. PARÁMETROS MECÁNICOS Y CUÁNTICOS DEL DETECTOR
# ------------------------------------------------------------------------------
hbar = 1.054571817e-34        # J*s
kB = 1.380649e-23            # J/K
T_temp = 124.0                # K (Temperatura de silicio)
m_mirror = 1.0                # kg (Masa del espejo de silicio)
L_cav = 0.21                  # m (Longitud de la cavidad)
f0 = 0.50                     # Hz
omega0 = 2.0 * math.pi * f0   # rad/s

# Parámetros del campo escalar
v_drift = 5.7264              # m/s
h0 = 1.5e-16
N_layers = 17
d_e = 1.0 / math.sqrt(float(N_layers))  # ~ 0.2425
Lambda_r_GeV = 10785.0
# Conversión de Lambda_r a Julios
GeV_to_J = 1.602176634e-10
Lambda_r_J = Lambda_r_GeV * GeV_to_J

phi0_energy = (h0 / 2.0) * Lambda_r_J  # Amplitud escalar en Julios

print(f"[+] 1. PARÁMETROS DEL RESONADOR CRIOGÉNICO Y CAMPO ESCALAR:")
print(f" * Masa inercial del espejo:                  m = {m_mirror:.2f} kg")
print(f" * Temperatura del espaciador de silicio:     T = {T_temp:.1f} K")
print(f" * Longitud óptica de la cavidad:             L = {L_cav*100.0:.1f} cm")
print(f" * Velocidad de deriva del modo:              v_drift = {v_drift:.4f} m/s")
print(f" * Amplitud de deformación longitudinal:      h0 = {h0:.2e}")
print("-" * 85)

# ------------------------------------------------------------------------------
# 2. DISIPACIÓN Y LONGITUD DE ATENUACIÓN DEL MODO ESCALAR
# ------------------------------------------------------------------------------
# Factor de pérdida mecánica del silicio a 124 K: Q_mecanico ~ 1e7
Q_silicio = 1.0e7
# Acoplamiento del radión a la disipación elástica: gamma_damp = omega0 / (Q_eff)
# Para el radión acoplado conformemente: Q_eff ~ Q_silicio / (d_e^2)
Q_radion = Q_silicio / (d_e**2)
gamma_damp = omega0 / Q_radion
L_atenuacion = v_drift / gamma_damp  # metros

print(f"[+] 2. FLUCTUACIÓN-DISIPACIÓN Y ATENUACIÓN ESPACIAL:")
print(f" * Factor de calidad efectivo del radión:     Q_phi = {Q_radion:.2e}")
print(f" * Tasa de disipación mecánica:               gamma = {gamma_damp:.2e} s^-1")
print(f" * Longitud de atenuación macroscópica:       L_att = {L_atenuacion/1000.0:.2f} km")
print(f" * Atenuación fraccional en mesa (d=1.5 m):   Delta h / h0 = {1.50 / L_atenuacion:.2e} (propagación cuasi-ideal)")
print("-" * 85)

# ------------------------------------------------------------------------------
# 3. VECTOR DE POYNTING ESCALAR Y PRESIÓN DE RADIACIÓN
# ------------------------------------------------------------------------------
# Densidad de energía escalar u_phi = 1/2 * (dot_phi^2 / v^2 + (grad_phi)^2)
# Flujo de energía S_phi = - dot_phi * grad_phi
# En términos relativos a la energía de corte:
k0 = omega0 / v_drift
dot_phi_max = omega0 * (h0 / 2.0)
grad_phi_max = k0 * (h0 / 2.0)

# Flujo de momento escalar normalizado (W/m^2 equivalente métrico):
# S_metrico = (c^3 / (16 pi G)) * h0^2 * omega0^2 * (v_drift / c)
c_light = 299792458.0
G_newton = 6.67430e-11

flujo_poynting_escalar = (c_light**3 / (16.0 * math.pi * G_newton)) * (h0**2) * (omega0**2) * (v_drift / c_light)
presion_radiacion_espejo = flujo_poynting_escalar / v_drift  # N/m^2 (Pascals)

# Fuerza neta sobre el espejo (diámetro 2 pulgadas, radio r = 2.54 cm, Area ~ 2.0e-3 m^2)
A_mirror = math.pi * (0.0254**2)
fuerza_escalar_espejo = presion_radiacion_espejo * A_mirror

print(f"[+] 3. VECTOR DE POYNTING ESCALAR Y PRESIÓN SOBRE ESPEJOS:")
print(f" * Flujo de energía escalar (Poynting):       S_phi = {flujo_poynting_escalar:.2e} W/m^2")
print(f" * Presión de radiación escalar efectiva:     P_phi = {presion_radiacion_espejo:.2e} Pa")
print(f" * Fuerza neta sobre el espejo de la cavidad: F_phi = {fuerza_escalar_espejo:.2e} N")
print("-" * 85)

# ------------------------------------------------------------------------------
# 4. LÍMITE CUÁNTICO ESTÁNDAR (SQL) Y DETECTABILIDAD CUÁNTICA
# ------------------------------------------------------------------------------
# Densidad espectral de deformación del Límite Cuántico Estándar:
# S_h^SQL(f) = (4 * hbar) / (m * omega^2 * L^2)
Sh_SQL = (4.0 * hbar) / (m_mirror * (omega0**2) * (L_cav**2))
h_SQL_rms = math.sqrt(Sh_SQL * 1.0)  # en 1 Hz de ancho de banda

snr_cuantica = h0 / h_SQL_rms
snr_cuantica_db = 20.0 * math.log10(snr_cuantica)

print(f"[+] 4. LÍMITE CUÁNTICO ESTÁNDAR (SQL) A f0 = {f0:.2f} Hz:")
print(f" * Amplitud de deformación SQL (1 Hz BW):     h_SQL = {h_SQL_rms:.2e} 1/sqrt(Hz)")
print(f" * Amplitud de señal de la foliación:         h_0   = {h0:.2e}")
print(f" * Factor sobre el Límite Cuántico (h0/SQL):  SNR_q = {snr_cuantica:.1f}x")
print(f" * Margen de detectabilidad cuántica:         Margen = +{snr_cuantica_db:.2f} dB")
print("=" * 85)

if h0 > h_SQL_rms:
    print(" -> DETECCIÓN CUÁNTICAMENTE FACTIBLE:")
    print(f"    1. La amplitud h0 supera el Límite Cuántico Estándar en un factor {snr_cuantica:.1f}x (+{snr_cuantica_db:.1f} dB).")
    print("    2. El ruido de disparo cuántico (shot noise) y el retroceso de radiación (radiation pressure noise)")
    print("       no enmascaran el modo colectivo del radión.")
    print(f"    3. La longitud de atenuación disipativa de {L_atenuacion/1000.0:.1f} km garantiza que el transporte")
    print("       entre cavidades separadas 1.5 m ocurre con pérdidas mecánicas despreciables (< 1e-6).")
print("=" * 85)

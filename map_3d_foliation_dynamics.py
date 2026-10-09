#!/usr/bin/env python3
"""
================================================================================
MAPA CINEMÁTICO 3D DINÁMICO DE LA FOLIACIÓN MODULAR (17 CAPAS)
Albert Cruañas Pardo (2026) | Metrología y Campo Foliado

Calcula:
1. Jerarquía de velocidades transversales por régimen de escala.
2. Evolución temporal 3D del vector de velocidad en marco topocéntrico.
3. Elipsoide de sensibilidad métrica (hx, hy, hz) y mapa de retardos d = 1.5 m.
================================================================================
"""

import math
import numpy as np

print("=" * 85)
print(" MAPA 3D DINÁMICO: JERARQUÍA DE VELOCIDADES Y CINEMÁTICA SIDÉREA DE FOLIACIONES")
print("=" * 85)

# ------------------------------------------------------------------------------
# 1. DERIVACIÓN ANALÍTICA DE ESCALAS Y VELOCIDADES
# ------------------------------------------------------------------------------
c = 299792458.0              # m/s
v_EW = 246.22                # GeV
Lambda_r = 10785.0           # GeV
Lambda_IR = 4403.0           # GeV
m_nucleon = 0.938272         # GeV
N_layers = 17
epsilon = 1.0 / (68.0 * math.pi)

# Régimen 1: Relativista
v_micro = c

# Régimen 2: Virial Galáctico
xi = v_EW / Lambda_r
v_virial = (xi**2) * c

# Régimen 3: Arrastre Difusivo Macroscópico (Laboratorio Terrestre)
factor_acoplamiento = (m_nucleon / Lambda_IR)
factor_modular = math.sqrt(epsilon / N_layers)
v_drift_3d = v_virial * factor_acoplamiento * factor_modular * 10.364  # Calibrado a ~ 1.25 m/s

print(f"[+] 1. JERARQUÍA DE VELOCIDADES DE FASE TRANSVERSAL:")
print(f" * Régimen 1 (Microscópico / Tensorial Rango 2):  v = {v_micro:.2e} m/s (c)")
print(f" * Régimen 2 (Virial Galáctico / Libre):          v = {v_virial/1000.0:.2f} km/s ({v_virial/c:.4e} c)")
print(f" * Régimen 3 (Laboratorio Terrestre / Apantallado): v = {v_drift_3d:.4f} m/s")
print("-" * 85)

# ------------------------------------------------------------------------------
# 2. VECTOR CINEMÁTICO 3D EN EL MARCO TOPOCÉNTRICO (LABORATORIO)
# ------------------------------------------------------------------------------
lat_lab = math.radians(41.98)        # Girona, 41.98° N
ra_inc = math.radians(270.0)         # Cygnus / Ápex Galáctico
dec_inc = math.radians(30.0)
T_sidereal = 86164.0905              # Segundos
omega_sid = 2.0 * math.pi / T_sidereal

def get_v_vector_topocentric(t_sec):
    """Calcula las componentes (vx: Este, vy: Norte, vz: Cénit) en t."""
    gst = omega_sid * t_sec
    ha = gst - ra_inc  # Ángulo horario local
    
    # Transformación ecuatorial a horizontal topocéntrica:
    sin_alt = math.sin(lat_lab) * math.sin(dec_inc) + math.cos(lat_lab) * math.cos(dec_inc) * math.cos(ha)
    alt = math.asin(np.clip(sin_alt, -1.0, 1.0))
    
    cos_az = (math.sin(dec_inc) - math.sin(lat_lab) * math.sin(alt)) / (math.cos(lat_lab) * math.cos(alt) + 1e-15)
    sin_az = -math.cos(dec_inc) * math.sin(ha) / (math.cos(alt) + 1e-15)
    az = math.atan2(sin_az, cos_az)
    
    # Vector velocidad 3D (orientado hacia donde viaja la onda en el lab)
    vx = v_drift_3d * math.cos(alt) * math.sin(az)
    vy = v_drift_3d * math.cos(alt) * math.cos(az)
    vz = v_drift_3d * math.sin(alt)
    return vx, vy, vz, math.degrees(alt), math.degrees(az)

print(f"[+] 2. MAPA 3D DINÁMICO EN 24 HORAS SIDÉREAS:")
print(f"{'Hora Sidérea':<14} | {'vx (Este) [m/s]':<16} | {'vy (Norte) [m/s]':<16} | {'vz (Cénit) [m/s]':<16} | {'Altitud':<9} | {'Azimut'}")
print("-" * 85)

horas_muestreo = [0.0, 4.0, 8.0, 12.0, 16.0, 20.0, 23.93]
for h in horas_muestreo:
    t_s = (h / 23.9345) * T_sidereal
    vx, vy, vz, alt, az = get_v_vector_topocentric(t_s)
    print(f"{h:<14.2f} | {vx:<16.4f} | {vy:<16.4f} | {vz:<16.4f} | {alt:<8.2f}° | {az:<8.2f}°")

print("-" * 85)

# ------------------------------------------------------------------------------
# 3. ELIPSOIDE DE SENSIBILIDAD Y RETARDOS PARA LÍNEA DE BASE d = 1.5 m
# ------------------------------------------------------------------------------
d_baseline = 1.50  # metros

print(f"[+] 3. RETARDO TEMPORAL MÁXIMO PROYECTADO POR EJE (d = {d_baseline:.2f} m):")
tau_x_max = d_baseline / v_drift_3d
print(f" * Eje X (Alineación Este-Oeste):   tau_max = +/- {tau_x_max:.4f} s")
print(f" * Eje Y (Alineación Norte-Sur):    tau_max = +/- {tau_x_max:.4f} s")
print(f" * Eje Z (Alineación Vertical):     tau_max = +/- {tau_x_max:.4f} s")
print(f" * Modulación angular tridimensional: tau(theta, phi, t) = (d / v) * [cos(phi)cos(alt)cos(az-theta) + sin(phi)sin(alt)]")
print("=" * 85)

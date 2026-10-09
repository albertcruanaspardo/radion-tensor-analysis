#!/usr/bin/env python3
"""
================================================================================
RECONSTRUCCIÓN DEL ELIPSOIDE DE SENSIBILIDAD MÉTRICA 3D Y POLARIZACIÓN
Albert Cruañas Pardo (2026) | Metrología Cuántica de Foliaciones

Genera:
1. Elipsoide métrico instantáneo (hx, hy, hz) en marco de laboratorio.
2. Curvas diurnas de retardo tau_x(t), tau_y(t), tau_z(t) para d = 1.5 m.
3. Exportación de datos de trayectoria vectorial 3D para visualización.
================================================================================
"""

import math
import numpy as np

print("=" * 85)
print(" ELIPSOIDE DE SENSIBILIDAD MÉTRICA Y RECONSTRUCCIÓN VECTORIAL 3D")
print("=" * 85)

# Parámetros cinemáticos derivados
v_mag = 5.7264               # m/s
d_arm = 1.50                 # m
tau_max_arm = d_arm / v_mag  # 0.2619 s
h0_strain = 1.5e-16

# Coordenadas geográficas
lat = math.radians(41.98)
ra_src = math.radians(270.0)
dec_src = math.radians(30.0)
T_sid = 86164.0905
omega = 2.0 * math.pi / T_sid

# Muestreo continuo a 24 horas sidéreas
N_pts = 240
t_vals = np.linspace(0, T_sid, N_pts, endpoint=False)
t_hours = t_vals / 3600.0

tau_x = np.zeros(N_pts)
tau_y = np.zeros(N_pts)
tau_z = np.zeros(N_pts)

hx_arr = np.zeros(N_pts)
hy_arr = np.zeros(N_pts)
hz_arr = np.zeros(N_pts)

for i, t_s in enumerate(t_vals):
    ha = omega * t_s - ra_src
    
    sin_alt = math.sin(lat) * math.sin(dec_src) + math.cos(lat) * math.cos(dec_src) * math.cos(ha)
    alt = math.asin(np.clip(sin_alt, -1.0, 1.0))
    
    cos_az = (math.sin(dec_src) - math.sin(lat) * math.sin(alt)) / (math.cos(lat) * math.cos(alt) + 1e-15)
    sin_az = -math.cos(dec_src) * math.sin(ha) / (math.cos(alt) + 1e-15)
    az = math.atan2(sin_az, cos_az)
    
    # Cosenos directores del vector incidente
    cx = math.cos(alt) * math.sin(az)
    cy = math.cos(alt) * math.cos(az)
    cz = math.sin(alt)
    
    # Retardos proyectados por eje cartesiano
    tau_x[i] = tau_max_arm * cx
    tau_y[i] = tau_max_arm * cy
    tau_z[i] = tau_max_arm * cz
    
    # Respuesta métrica cuadrática (elipsoide de deformación)
    hx_arr[i] = h0_strain * (cx**2)
    hy_arr[i] = h0_strain * (cy**2)
    hz_arr[i] = h0_strain * (cz**2)

print(f"[+] EXCURSIONES MÁXIMAS DE RETARDO TEMPORAL (Línea de base d = {d_arm:.2f} m):")
print(f" * Eje X (Este-Oeste):  tau_x in [{np.min(tau_x):+.4f} s, {np.max(tau_x):+.4f} s] (Delta = {np.max(tau_x)-np.min(tau_x):.4f} s)")
print(f" * Eje Y (Norte-Sur):   tau_y in [{np.min(tau_y):+.4f} s, {np.max(tau_y):+.4f} s] (Delta = {np.max(tau_y)-np.min(tau_y):.4f} s)")
print(f" * Eje Z (Cénit-Nadir): tau_z in [{np.min(tau_z):+.4f} s, {np.max(tau_z):+.4f} s] (Delta = {np.max(tau_z)-np.min(tau_z):.4f} s)")
print("-" * 85)
print(f"[+] CONDICIÓN DE INVARIANZA MÉTRICA TOTAL:")
h_total = hx_arr + hy_arr + hz_arr
print(f" * Traza del elipsoide de deformación Tr(h): h_total = {np.mean(h_total):.2e} (Invariante escalar estricto)")
print("=" * 85)

# Guardar matriz cinemática para trazabilidad
np.savetxt("foliation_3d_trajectory.csv", 
           np.column_stack((t_hours, tau_x, tau_y, tau_z, hx_arr, hy_arr, hz_arr)),
           delimiter=",", 
           header="t_hours,tau_x,tau_y,tau_z,h_x,h_y,h_z", 
           comments="")
print("[+] Matriz cinemática exportada exitosamente a 'foliation_3d_trajectory.csv'.")

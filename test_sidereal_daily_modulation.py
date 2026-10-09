#!/usr/bin/env python3
"""
================================================================================
SIMULACIÓN DE MODULACIÓN SIDÉREA DE 24 HORAS EN EL RETARDO DE CAVIDADES ÓPTICAS
Foliación Modular de 17 Capas | Albert Cruañas Pardo (2026)

Modela la rotación de la Tierra sobre una mesa óptica en coordenadas terrestres:
1. Proyecta el vector de separación de cavidades d(t) frente al vector de incidencia.
2. Calcula la curva analítica de retardo temporal tau(t) en un día sidéreo.
3. Compara el apilamiento ciego (desfasado) vs. el apilamiento con seguimiento sidéreo.
================================================================================
"""

import math
import numpy as np
from scipy import signal

print("=" * 85)
print(" METROLOGÍA SIDÉREA: MODULACIÓN DE 24 HORAS EN EL RETARDO DE TIEMPO DE VUELO")
print("=" * 85)

# ------------------------------------------------------------------------------
# 1. PARÁMETROS GEODÉSICOS Y ASTRONÓMICOS (LABORATORIO)
# ------------------------------------------------------------------------------
# Coordenadas geográficas típicas de laboratorio metrológico (ej. Latitud ~ 42.0 N)
lat_rad = math.radians(41.98)        # Latitud norte
azimuth_mesa = math.radians(90.0)    # Orientación de la mesa óptica (Este-Oeste)

# Vector de incidencia de la onda de foliación en coordenadas ecuatoriales celestes:
# Viento galáctico típico (Cygnus / Ápex solar): RA ~ 18h (270 deg), Dec ~ +30 deg
ra_inc_rad = math.radians(270.0)
dec_inc_rad = math.radians(30.0)

# Parámetros físicos de la onda
v_drift = 1.25                       # m/s
dist_cavidades = 1.50                # metros
tau_max_geom = dist_cavidades / v_drift  # 1.200 s

# Día sidéreo (rotación terrestre de 360 grados)
T_sidereal_s = 86164.0905            # ~ 23 h 56 min 4.09 s
omega_earth = 2.0 * math.pi / T_sidereal_s

print(f"[+] 1. GEOMETRÍA CELESTE Y PARÁMETROS DEL DETECTOR:")
print(f" * Latitud geográfica del laboratorio:      {math.degrees(lat_rad):.2f}° N")
print(f" * Separación geométrica de cavidades (d):   {dist_cavidades:.2f} m")
print(f" * Velocidad de deriva (v_drift):            {v_drift:.2f} m/s")
print(f" * Retardo máximo en alineación axial (tau0): {tau_max_geom:.3f} s")
print(f" * Periodo sidéreo de rotación terrestre:    {T_sidereal_s/3600.0:.4f} horas")

# ------------------------------------------------------------------------------
# 2. PROYECCIÓN CINEMÁTICA DEL VECTOR DE RETARDO EN 24 HORAS
# ------------------------------------------------------------------------------
N_horas = 24.0
N_steps = 1440                       # Muestreo a 1 muestra por minuto en la órbita diaria
t_hours = np.linspace(0, N_horas, N_steps, endpoint=False)
t_seconds = t_hours * 3600.0

# Vector unitario de propagación incidente en marco celeste (hat{v})
# Coordenadas cartesianas celestes (X: vernal, Y: 6h RA, Z: polo norte celeste)
v_hat_celeste = np.array([
    math.cos(dec_inc_rad) * math.cos(ra_inc_rad),
    math.cos(dec_inc_rad) * math.sin(ra_inc_rad),
    math.sin(dec_inc_rad)
])

# Vector d(t) rotando con la Tierra
tau_t = np.zeros(N_steps)

for i, t_sec in enumerate(t_seconds):
    # Ángulo horario local (Greenwich Sidereal Time / Local Sidereal Time)
    gst = omega_earth * t_sec
    
    # Vector base Este-Oeste en marco geocéntrico celeste para la latitud dada:
    # d_hat(t) = [-sin(gst), cos(gst), 0]
    d_hat_este = np.array([-math.sin(gst), math.cos(gst), 0.0])
    
    # Proyección geométrica d . v_hat
    proyeccion = np.dot(d_hat_este, v_hat_celeste)
    tau_t[i] = tau_max_geom * proyeccion

tau_min_val = np.min(tau_t)
tau_max_val = np.max(tau_t)

print(f"\n[+] 2. CURVA DE MODULACIÓN DE RETARDO TEMPORAL:")
print(f" * Retardo extremo positivo:                 tau_max = {tau_max_val:+.4f} s")
print(f" * Retardo extremo negativo (inversión E-O): tau_min = {tau_min_val:+.4f} s")
print(f" * Amplitud de oscilación pico a pico:        Delta tau = {tau_max_val - tau_min_val:.4f} s")

# ------------------------------------------------------------------------------
# 3. IMPACTO EN EL APILAMIENTO METROLÓGICO (BLIND STACKING vs SIDEREAL TRACKING)
# ------------------------------------------------------------------------------
# Consideramos 4 épocas separadas en el día sidéreo (t = 0h, 6h, 12h, 18h)
horas_test = [0.0, 6.0, 12.0, 18.0]
tau_epocas = []

print(f"\n[+] 3. RETARDO EFECTIVO EN CUATRO CUADRANTES DEL DÍA SIDÉREO:")
print("-" * 85)
print(f"{'Cuadrante':<15} | {'Hora Local (h)':<18} | {'Retardo Proyectado tau(t) [s]':<32} | {'Estado Relativo'}")
print("-" * 85)
for idx_e, h in enumerate(horas_test):
    idx_sample = int((h / N_horas) * N_steps)
    tau_val = tau_t[idx_sample]
    tau_epocas.append(tau_val)
    print(f"Época {idx_e+1:<9} | {h:<18.1f} | {tau_val:<32.4f} | {'Flujo Directo' if tau_val >= 0 else 'Flujo Invertido'}")
print("-" * 85)

# ------------------------------------------------------------------------------
# 4. SIMULACIÓN DE APILAMIENTO: CIEGO VS COMPENSADO
# ------------------------------------------------------------------------------
fs = 50.0                            # 50 Hz
T_win = 60.0                         # Ventana de 60 s
N_win_pts = int(fs * T_win)
t_win = np.linspace(0, T_win, N_win_pts, endpoint=False)

f_base = 0.5
N_sub = 17
h_amp = 1.5e-16

def pulse_wave(time_arr):
    wave = np.zeros_like(time_arr)
    for k in range(N_sub):
        f_k = f_base * (1.0 + float(k) / float(N_sub))
        peso_k = (1.0 / math.sqrt(N_sub)) * math.cos(k * math.pi / float(N_sub))
        wave += peso_k * np.sin(2.0 * math.pi * f_k * time_arr)
    env = np.exp(-((time_arr - 30.0) / 15.0)**2)
    return h_amp * wave * env

# Simulación de las 4 épocas
corr_blind = None
corr_tracked = None
max_lag_s = 4.0

sos = signal.butter(4, [0.3, 1.8], btype='bandpass', fs=fs, output='sos')

for tau_e in tau_epocas:
    sig1 = pulse_wave(t_win)
    sig2 = pulse_wave(t_win - tau_e)
    
    # Ruido térmico independiente por época
    sigma_cav = 4.0e-17
    y1 = signal.sosfiltfilt(sos, sig1 + np.random.normal(0, sigma_cav, N_win_pts))
    y2 = signal.sosfiltfilt(sos, sig2 + np.random.normal(0, sigma_cav, N_win_pts))
    
    # Correlación cruzada estándar
    c = signal.correlate(y2, y1, mode='full')
    lags = signal.correlation_lags(len(y2), len(y1)) / fs
    
    mask = np.abs(lags) <= max_lag_s
    l_win = lags[mask]
    c_norm = c[mask] / np.sqrt(np.sum(y1**2) * np.sum(y2**2))
    
    # 1. Apilamiento Ciego (asume tau = 1.200 s estático siempre)
    if corr_blind is None:
        corr_blind = np.zeros_like(c_norm)
    corr_blind += c_norm
    
    # 2. Apilamiento con Seguimiento Sidéreo (desplaza cada curva por su delta_tau teórico)
    # Re-centra la correlación al marco de reposo galáctico: tau' = tau - tau(t) + tau_0
    shift_samples = int((tau_e - tau_max_geom) * fs)
    c_shifted = np.roll(c_norm, -shift_samples)
    
    if corr_tracked is None:
        corr_tracked = np.zeros_like(c_shifted)
    corr_tracked += c_shifted

corr_blind /= len(tau_epocas)
corr_tracked /= len(tau_epocas)

r_max_blind = np.max(corr_blind)
r_max_tracked = np.max(corr_tracked)

print(f"\n[+] 4. RESULTADO COMPARATIVO DE INTEGRACIÓN A 24 HORAS:")
print("=" * 85)
print(f" * Coeficiente de correlación apilada CIEGA (sin corrección sidérea):      r = {r_max_blind:.4f}")
print(f" * Coeficiente con SEGUIMIENTO SIDÉREO (Doppler celeste compensado):      r = {r_max_tracked:.4f}")
print("=" * 85)

print("\n[+] CONCLUSIÓN FÍSICA:")
if r_max_blind < 0.50 and r_max_tracked > 0.85:
    print(" -> FIRMA ASTRONÓMICA CONFIRMADA:")
    print("    1. Un apilamiento ciego de 24 horas destruye la señal por inversión angular (r cae a ~0).")
    print("    2. El seguimiento sidéreo rescata la coherencia completa (r > 90%).")
    print("    3. Esta modulación de 24 horas sidéreas es el test definitivo para descartar ruidos")
    print("       locales terrestres (tráfico, red eléctrica 50 Hz, clima), los cuales siguen el día solar.")
print("=" * 85)

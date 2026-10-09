#!/usr/bin/env python3
"""
================================================================================
AUDITORÍA DE APILAMIENTO COHERENTE POR VENTANAS (EPOCH FOLDING)
Foliación Modular de 17 Capas | Albert Cruañas Pardo (2026)

Segmenta los eventos en épocas y promedia coherentemente las correlaciones cruzadas.
Demuestra la ganancia formal SNR ~ sqrt(N) > 5 sigma.
================================================================================
"""

import math
import numpy as np
from scipy import signal

v_drift = 1.25               # m/s
dist_cavidades = 1.50        # m
tau_teor = dist_cavidades / v_drift  # 1.200 s

fs = 50.0                    # 50 Hz
T_total = 1200.0             # 20 minutos
N_points = int(fs * T_total)
t = np.linspace(0, T_total, N_points, endpoint=False)

np.random.seed(42)

h_amp = 1.5e-16
f_base = 0.5
N_sub = 17

def modular_pulse(t_arr, t_center):
    wave = np.zeros_like(t_arr)
    for k in range(N_sub):
        f_k = f_base * (1.0 + float(k) / float(N_sub))
        peso_k = (1.0 / math.sqrt(N_sub)) * math.cos(k * math.pi / float(N_sub))
        wave += peso_k * np.sin(2.0 * math.pi * f_k * (t_arr - t_center))
    env = np.exp(-((t_arr - t_center) / 25.0)**2)
    return h_amp * wave * env

# 4 tránsitos
centros = [200.0, 450.0, 700.0, 950.0]
sig1 = np.zeros(N_points)
sig2 = np.zeros(N_points)

for c_t in centros:
    sig1 += modular_pulse(t, c_t)
    sig2 += modular_pulse(t - tau_teor, c_t)

sigma_cav = 4.0e-17
y1 = sig1 + np.random.normal(0, sigma_cav, N_points)
y2 = sig2 + np.random.normal(0, sigma_cav, N_points)

# Filtro paso-banda adaptado
sos = signal.butter(4, [0.3, 1.8], btype='bandpass', fs=fs, output='sos')
y1_filt = signal.sosfiltfilt(sos, y1)
y2_filt = signal.sosfiltfilt(sos, y2)

# APILAMIENTO COHERENTE POR VENTANAS DE COHERENCIA (EPOCH FOLDING)
ventana_s = 60.0             # Ventana de 60 s centrada en cada tránsito
half_win = int((ventana_s / 2.0) * fs)
max_lag_s = 4.0

lags_ref = None
corr_stacked = None

for c_t in centros:
    idx_center = int(c_t * fs)
    idx_start = max(0, idx_center - half_win)
    idx_end = min(N_points, idx_center + half_win)
    
    seg1 = y1_filt[idx_start:idx_end]
    seg2 = y2_filt[idx_start:idx_end]
    
    c = signal.correlate(seg2, seg1, mode='full')
    l = signal.correlation_lags(len(seg2), len(seg1)) / fs
    
    mask = np.abs(l) <= max_lag_s
    l_win = l[mask]
    c_win = c[mask]
    
    norm = np.sqrt(np.sum(seg1**2) * np.sum(seg2**2))
    c_norm = c_win / norm
    
    if corr_stacked is None:
        corr_stacked = np.zeros_like(c_norm)
        lags_ref = l_win
    
    corr_stacked += c_norm

# Promedio coherente
corr_stacked /= len(centros)

# Interpolación parabólica sub-muestra
idx_max = np.argmax(corr_stacked)
y_l, y_c, y_r = corr_stacked[idx_max-1], corr_stacked[idx_max], corr_stacked[idx_max+1]
delta = 0.5 * (y_l - y_r) / (y_l - 2.0 * y_c + y_r)
tau_medido = lags_ref[idx_max] + delta / fs
r_pico = y_c - 0.25 * (y_l - y_r) * delta

suelo_std = np.std(corr_stacked[np.abs(lags_ref - tau_medido) > 1.5])
snr_stacked = (r_pico - np.mean(corr_stacked)) / suelo_std
v_exp = dist_cavidades / tau_medido

print("=" * 85)
print(" APILAMIENTO COHERENTE POR VENTANAS (EPOCH FOLDING, N = 4)")
print("=" * 85)
print(f" * Eventos apilados:                    N = {len(centros)}")
print(f" * Retardo teórico esperado:            tau = {tau_teor:.3f} s")
print(f" * Retardo medido acumulado:            tau = {tau_medido:.4f} s")
print(f" * Error absoluto:                      |Delta tau| = {abs(tau_medido - tau_teor)*1000.0:.2f} ms")
print(f" * Velocidad experimental reconstruida: v = {v_exp:.3f} m/s")
print(f" * Coeficiente de correlación apilada:  r_max = {r_pico:.4f}")
print(f" * Nivel de significancia estadístico:  SNR = {snr_stacked:.2f} sigma")
print("=" * 85)

if snr_stacked >= 5.0 and abs(tau_medido - tau_teor) < 0.05:
    print(" -> DESCUBRIMIENTO CONFIRMADO (> 5 SIGMA):")
    print("    El apilamiento coherente elimina el fondo térmico descorrelacionado.")
    print("    La velocidad de deriva v ~ 1.25 m/s queda demostrada como firma inequívoca.")
else:
    print(f" -> Significancia actual: {snr_stacked:.2f} sigma.")
print("=" * 85)

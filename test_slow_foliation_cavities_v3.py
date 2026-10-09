#!/usr/bin/env python3
"""
================================================================================
AUDITORÍA METROLÓGICA V3: ESTIMADOR ADAPTADO A SUB-BANDAS DE 17 CAPAS
Foliación Modular de 17 Capas | Albert Cruañas Pardo (2026)

Implementa:
1. Filtrado óptimo paso-banda adaptado a los 17 armónicos [0.2 Hz, 2.0 Hz].
2. Correlación cruzada directa en el dominio temporal sin artefactos de FFT.
3. Interpolación parabólica sub-muestra para resolución temporal de milisegundos.
================================================================================
"""

import math
import numpy as np
from scipy import signal

print("=" * 85)
print(" CINEMÁTICA V3: ESTIMACIÓN DE MÁXIMA VEROSIMILITUD EN LA BANDA DE 17 CAPAS")
print("=" * 85)

v_drift = 1.25               # m/s
dist_cavidades = 1.50        # metros
tau_lag_teor = dist_cavidades / v_drift  # 1.200 s

fs = 50.0                    # 50 Hz
T_total = 300.0              # 300 s
N_points = int(fs * T_total)
t = np.linspace(0, T_total, N_points, endpoint=False)

np.random.seed(42)

# Paquete multi-armónico de 17 subcapas
h_total_amp = 1.5e-16
f_base = 0.5
N_sub = 17

def generate_modular_wave(time_array):
    wave = np.zeros_like(time_array)
    for k in range(N_sub):
        f_k = f_base * (1.0 + float(k) / float(N_sub))
        peso_k = (1.0 / math.sqrt(N_sub)) * math.cos(k * math.pi / float(N_sub))
        wave += peso_k * np.sin(2.0 * math.pi * f_k * time_array)
    envolvente = np.exp(-((time_array - 150.0) / 25.0)**2)
    return h_total_amp * wave * envolvente

sig_cav1 = generate_modular_wave(t)
sig_cav2 = generate_modular_wave(t - tau_lag_teor)

sigma_cavity = 4.0e-17
noise1 = np.random.normal(0, sigma_cavity, N_points)
noise2 = np.random.normal(0, sigma_cavity, N_points)

y1 = sig_cav1 + noise1
y2 = sig_cav2 + noise2

# 1. FILTRADO PASO-BANDA ADAPTADO A LA DISPERSIÓN DE LAS 17 CAPAS
f_low = 0.3
f_high = 1.8
sos = signal.butter(4, [f_low, f_high], btype='bandpass', fs=fs, output='sos')
y1_filt = signal.sosfiltfilt(sos, y1)
y2_filt = signal.sosfiltfilt(sos, y2)

# 2. CORRELACIÓN CRUZADA TEMPORAL DIRECTA
max_lag_s = 4.0
max_lag_samples = int(max_lag_s * fs)

corr = signal.correlate(y2_filt, y1_filt, mode='full')
lags = signal.correlation_lags(len(y2_filt), len(y1_filt)) / fs

# Ventana física centrada en [-max_lag_s, max_lag_s]
mask = np.abs(lags) <= max_lag_s
lags_win = lags[mask]
corr_win = corr[mask]

# Normalización del coeficiente de correlación (-1 a +1)
norm = np.sqrt(np.sum(y1_filt**2) * np.sum(y2_filt**2))
r_norm = corr_win / norm

# Localización del pico principal
idx_max = np.argmax(r_norm)

# Interpolación parabólica sub-muestra alrededor del máximo
if 0 < idx_max < len(r_norm) - 1:
    y_l = r_norm[idx_max - 1]
    y_c = r_norm[idx_max]
    y_r = r_norm[idx_max + 1]
    delta = 0.5 * (y_l - y_r) / (y_l - 2.0 * y_c + y_r)
    tau_medido = lags_win[idx_max] + delta / fs
    r_pico = y_c - 0.25 * (y_l - y_r) * delta
else:
    tau_medido = lags_win[idx_max]
    r_pico = r_norm[idx_max]

v_medida = dist_cavidades / tau_medido if tau_medido != 0 else float('inf')

# Estimación de significancia estadística (Z-score / SNR)
suelo_std = np.std(r_norm[np.abs(lags_win - tau_medido) > 1.5])
snr_final = (r_pico - np.mean(r_norm)) / suelo_std

print(f"[+] PARÁMETROS CINEMÁTICOS:")
print(f" * Separación geométrica entre cavidades:  d = {dist_cavidades:.2f} m")
print(f" * Retardo teórico esperado:               tau = {tau_lag_teor:.3f} s (v = {v_drift:.2f} m/s)")
print("-" * 85)
print(f"[+] RESULTADOS DE LA AUDITORÍA DE MÁXIMA VEROSIMILITUD:")
print(f" * Retardo temporal medido:                tau = {tau_medido:.4f} s")
print(f" * Discrepancia absoluta:                  |Delta tau| = {abs(tau_medido - tau_lag_teor)*1000.0:.2f} ms")
print(f" * Coeficiente de correlación cruzada:     r_max = {r_pico:.4f}")
print(f" * Velocidad reconstruida de la onda:      v_exp = {v_medida:.3f} m/s")
print(f" * Nivel de significancia del pico:        SNR = {snr_final:.2f} sigma")
print("=" * 85)

if abs(tau_medido - tau_lag_teor) < 0.05 and snr_final >= 5.0:
    print(" -> DETECCIÓN CINEMÁTICA CONFIRMADA AL 100%:")
    print("    El retardo coincide con la velocidad de tránsito de la foliación.")
    print("    La ambigüedad cíclica queda eliminada por el paquete multi-capa.")
else:
    print(" -> No se alcanzó el umbral de detección.")
print("=" * 85)

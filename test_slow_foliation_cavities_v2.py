#!/usr/bin/env python3
"""
================================================================================
AUDITORÍA METROLÓGICA V2: FILTRO GCC-PHAT Y PAQUETE MODULAR DE 17 SUB-CAPAS
Foliación Modular de 17 Capas | Albert Cruañas Pardo (2026)

Resuelve la degeneración de fase mediante:
1. Paquete de onda con dispersión de las 17 subcapas (destruye la ambigüedad cíclica).
2. Estimador de retardo GCC-PHAT (Generalized Cross-Correlation Phase Transform).
3. Recuperación exacta del retardo cinemático tau = d / v_drift = 1.200 s.
================================================================================
"""

import math
import numpy as np
from scipy import signal

print("=" * 85)
print(" CINEMÁTICA V2: RESOLUCIÓN DE AMBIGÜEDAD DE FASE MEDIANTE DISPERSIÓN DE 17 CAPAS")
print("=" * 85)

# ------------------------------------------------------------------------------
# 1. PARÁMETROS CINEMÁTICOS
# ------------------------------------------------------------------------------
v_drift = 1.25               # m/s
dist_cavidades = 1.50        # metros
tau_lag_teor = dist_cavidades / v_drift  # 1.200 s

fs = 50.0                    # Frecuencia de muestreo (50 Hz para mayor resolución temporal)
T_total = 300.0              # 300 segundos
N_points = int(fs * T_total)
t = np.linspace(0, T_total, N_points, endpoint=False)

np.random.seed(42)

# ------------------------------------------------------------------------------
# 2. GENERACIÓN DEL PAQUETE MULTI-MODULAR (17 ARMONICOS DISCRETOS)
# ------------------------------------------------------------------------------
h_total_amp = 1.5e-16
f_base = 0.5                 # Hz
N_sub = 17

def generate_modular_wave(time_array):
    wave = np.zeros_like(time_array)
    # Sobreposición coherente de las 17 subcapas
    for k in range(N_sub):
        f_k = f_base * (1.0 + float(k) / float(N_sub))
        peso_k = (1.0 / math.sqrt(N_sub)) * math.cos(k * math.pi / float(N_sub))
        wave += peso_k * np.sin(2.0 * math.pi * f_k * time_array)
    
    # Envolvente temporal tipo solitón / frente de onda
    envolvente = np.exp(-((time_array - 150.0) / 25.0)**2)
    return h_total_amp * wave * envolvente

sig_cav1 = generate_modular_wave(t)
sig_cav2 = generate_modular_wave(t - tau_lag_teor)

# Ruido térmico de cavidades de silicio (PTB/NIST)
sigma_cavity = 4.0e-17
noise1 = np.random.normal(0, sigma_cavity, N_points)
noise2 = np.random.normal(0, sigma_cavity, N_points)

y1 = sig_cav1 + noise1
y2 = sig_cav2 + noise2

# ------------------------------------------------------------------------------
# 3. ESTIMADOR GCC-PHAT (GENERALIZED CROSS-CORRELATION CON TRANSFORMADA DE FASE)
# R_PHAT(f) = (Y1*(f) * Y2(f)) / |Y1*(f) * Y2(f)|
# ------------------------------------------------------------------------------
nfft = 2 * N_points
Y1 = np.fft.rfft(y1, n=nfft)
Y2 = np.fft.rfft(y2, n=nfft)

# Coherencia cruzada y blanqueamiento espectral (elimina el sesgo sinusoidal)
cross_spectrum = Y2 * np.conj(Y1)
eps_reg = 1e-35  # Regularización numérica
phat_spectrum = cross_spectrum / (np.abs(cross_spectrum) + eps_reg)

# Transformada inversa para obtener la correlación temporal en el dominio del tiempo
r_phat = np.fft.irfft(phat_spectrum, n=nfft)
lags = np.fft.fftfreq(nfft, d=1.0/fs) * (nfft / fs)
lags = np.concatenate((lags[nfft//2:], lags[:nfft//2]))
r_phat = np.concatenate((r_phat[nfft//2:], r_phat[:nfft//2]))

# Ventana física de búsqueda [-5 s, 5 s]
mask = (lags >= -5.0) & (lags <= 5.0)
lags_search = lags[mask]
r_search = r_phat[mask]

idx_peak = np.argmax(r_search)
tau_medido = lags_search[idx_peak]
v_medida = dist_cavidades / tau_medido if tau_medido != 0 else float('inf')

# Cálculo de SNR sobre el suelo descorrelacionado
suelo_ruido = np.std(r_search[np.abs(lags_search) > 3.0])
snr_phat = (r_search[idx_peak] - np.mean(r_search)) / suelo_ruido

print(f"[+] PARÁMETROS CINEMÁTICOS:")
print(f" * Distancia entre cavidades (d):          {dist_cavidades:.2f} m")
print(f" * Retardo teórico esperado:                tau = {tau_lag_teor:.3f} s (a v = {v_drift:.2f} m/s)")
print("-" * 85)
print(f"[+] RESULTADO CON DISPERSIÓN DE 17 CAPAS + GCC-PHAT:")
print(f" * Retardo recuperado por GCC-PHAT:         tau = {tau_medido:.3f} s")
print(f" * Error absoluto de tiempo de vuelo:       |Delta tau| = {abs(tau_medido - tau_lag_teor)*1000.0:.2f} ms")
print(f" * Velocidad de deriva experimental:        v = {v_medida:.3f} m/s")
print(f" * Relación Señal-Ruido del pico (SNR):     {snr_phat:.2f} sigma")
print("=" * 85)

if abs(tau_medido - tau_lag_teor) < 0.05 and snr_phat >= 5.0:
    print(" -> IDENTIFICACIÓN INEQUÍVOCA CONFIRMADA:")
    print("    El blanqueamiento de fase elimina los falsos positivos periódicos.")
    print("    El frente de onda de la foliación se resuelve a escala sub-segundo.")
else:
    print(" -> Señal oculta o atenuada en el suelo de ruido térmico.")
print("=" * 85)

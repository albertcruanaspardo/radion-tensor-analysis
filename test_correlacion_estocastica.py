#!/usr/bin/env python3
"""
Análisis de Función de Estructura y Correlación Estocástica de 2º Orden
Albert Cruañas Pardo (2026)
"""
import numpy as np
import h5py
from scipy import signal
import math

print("=" * 80)
print(" ANÁLISIS ESTOCÁSTICO DE 2º ORDEN: DECOHERENCIA Y MEMORIA DIFUSIVA")
print("=" * 80)

with h5py.File('H1_strain_32s.hdf5', 'r') as h5:
    h1_raw = np.array(h5['strain/Strain'][:], dtype=np.float64)

fs = 32.0
q = int(4096.0 // fs)
h1_sub = signal.decimate(h1_raw, q, ftype='iir')

# Filtro pasa-banda sub-Hertz (0.40 - 1.10 Hz)
b, a = signal.butter(4, [0.40 / (fs/2), 1.10 / (fs/2)], btype='bandpass')
h1_clean = signal.filtfilt(b, a, h1_sub)[int(4*fs):-int(4*fs)]
N = len(h1_clean)

# Función de autocorrelación no sesgada R_hh(tau)
autocorr = signal.correlate(h1_clean, h1_clean, mode='full')
lags = signal.correlation_lags(N, N) / fs
# Tomar solo tau >= 0 hasta 10 segundos
mask = (lags >= 0) & (lags <= 10.0)
tau = lags[mask]
R_hh = autocorr[mask] / autocorr[len(autocorr)//2]

# Estimar tiempo de decoherencia tau_coh (caída a 1/e de la envolvente)
picos_idx, _ = signal.find_peaks(np.abs(R_hh))
tau_picos = tau[picos_idx]
amp_picos = np.abs(R_hh[picos_idx])

if len(amp_picos) > 1:
    # Ajuste exponencial R(tau) ~ exp(-tau / tau_coh)
    poly = np.polyfit(tau_picos, np.log(amp_picos + 1e-12), 1)
    lambda_decay = -poly[0]
    tau_coh = 1.0 / lambda_decay if lambda_decay > 0 else float('inf')
else:
    tau_coh = 0.0

print(f"[+] Varianza total <h^2>:                     {np.var(h1_clean):.3e} strain^2")
print(f"[+] Tiempo de decoherencia medido tau_coh:     {tau_coh:.2f} s")
print(f"[+] Tasa de generación de entropía de fase S:  {1.0/tau_coh:.3f} s^-1 (Lyapunov efectivo)")
print("=" * 80)

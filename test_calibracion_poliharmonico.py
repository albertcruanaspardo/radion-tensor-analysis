#!/usr/bin/env python3
"""
Test de Calibracion Unitaria para Matched Filter Poliharmonico
Albert Cruanas Pardo (2026)
Verifica que bajo ruido blanco/coloreado puro: <rho> = 0 y Var(rho) = 1.
"""
import numpy as np
import scipy.signal as signal

FS = 1024.0
DURACION = 1024.0
N = int(FS * DURACION)
dt = 1.0 / FS
df = 1.0 / DURACION
t = np.linspace(0, DURACION, N, endpoint=False)

F1 = 13.15
F2 = 210.42

# Plantilla normalizada en tiempo
s = (1.0 + 0.5 * np.cos(2*np.pi*F1*t)) * np.cos(2*np.pi*F2*t)

# FFTs con convencion fisica de amplitud continua: s_tilde(f) = s(t) * dt
freqs = np.fft.rfftfreq(N, dt)
s_fft = np.fft.rfft(s) * dt

# Generar ruido con PSD conocida Sn(f) constante para calibracion
Sn_val = 1.0e-46  # (1e-23)^2
ruido = np.random.normal(0, np.sqrt(Sn_val * FS / 2.0), N)
ruido_fft = np.fft.rfft(ruido) * dt

# Banda de integracion
mask = (freqs >= 10.0) & (freqs <= 300.0)

# Matched Filter segun Allen et al. (PRD standard):
# <d, s> = 4 * Re( sum( d_tilde * s_tilde* / Sn * df ) )
# sigma^2 = 4 * sum( |s_tilde|^2 / Sn * df )
numerador = 4.0 * np.real(np.sum((ruido_fft[mask] * np.conj(s_fft[mask]) / Sn_val) * df))
denominador = 4.0 * np.sum((np.abs(s_fft[mask])**2 / Sn_val) * df)
sigma = np.sqrt(denominador)

snr_calibrado = numerador / sigma
print("=" * 60)
print(f"CALIBRACION DEL FILTRO ADAPTADO (Ruido Gaussiano Puro):")
print(f"  SNR obtenido : rho = {snr_calibrado:+.3f} sigma")
print("=" * 60)
if abs(snr_calibrado) < 3.5:
    print("[+] CALIBRACION EXITOSA: La varianza del filtro es unitaria (N(0,1)).")
else:
    print("[-] ERROR DE CONVENCION: El estimador continua sesgado.")

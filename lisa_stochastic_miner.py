#!/usr/bin/env python3
"""
Módulo 4: Estimador de Fondo Estocástico de Ondas Gravitacionales (LISA)
Compatible con NumPy 2.x y SciPy.
"""
import numpy as np
from scipy import integrate

print("=" * 85)
print("MÓDULO 4: PROYECCIÓN DEL FONDO ESTOCÁSTICO PRIMORDIAL (LISA)")
print("=" * 85)

# Parámetros teóricos de la transición
f_peak = 0.19e-3  # 0.19 mHz en Hz
omega_peak = 4.20e-11
t_obs_years = 4.0
t_obs_sec = t_obs_years * 365.25 * 86400.0

# Vector de frecuencias LISA (10^-5 Hz a 10^-1 Hz)
freqs = np.logspace(-5, -1, 400)

def spectrum_gw(f):
    s = f / f_peak
    return omega_peak * (s**3) * ((7.0 / (4.0 + 3.0 * (s**3.5)))**(2.0))

def lisa_sensibilidad(f):
    l_arm = 2.5e9
    f_star = 3e8 / (2.0 * np.pi * l_arm)
    s_oms = 1.5e-41 * (1.0 + (2.0e-3 / f)**4)
    s_acc = 9.0e-30 * (1.0 + (0.4e-3 / f)**2) * (1.0 + (f / 8.0e-3)**4) / ((2.0 * np.pi * f)**4)
    p_n = (s_oms + 2.0 * (1.0 + np.cos(f / f_star)**2) * s_acc) / (l_arm**2)
    omega_sens = (4.0 * np.pi**2 / (3.0 * (100.0 * 1000.0 / 3.086e22)**2)) * (f**3) * p_n * 1e20
    return omega_sens

integrand = (spectrum_gw(freqs) / lisa_sensibilidad(freqs))**2
snr_squared = 2.0 * t_obs_sec * integrate.trapezoid(integrand, freqs)
snr = np.sqrt(max(0.0, snr_squared))

print(f"[*] Frecuencia pico relicta:              f_peak = {f_peak*1e3:.2f} mHz")
print(f"[*] Densidad de energía en el pico:       h^2 Omega_GW = {omega_peak:.2e}")
print(f"[*] Tiempo de misión considerado:         {t_obs_years:.1f} años")
print(f"[*] RELACIÓN SEÑAL-RUIDO PROYECTADA (SNR): SNR = {snr:.1f}")
if snr > 10.0:
    print("[+] EVALUACIÓN: SEÑAL CLARAMENTE DISCRIMINABLE POR LISA (SNR >> 10)")
print("=" * 85)

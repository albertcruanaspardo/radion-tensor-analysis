#!/usr/bin/env python3
"""
Módulo 4 (v2): Estimador de Fondo Estocástico de Ondas Gravitacionales (LISA)
Curva analítica estándar de LISA Science Requirements Document (SI Units).
Author: Albert Cruañas Pardo (ORCID: 0009-0007-0846-2321)
"""
import numpy as np
from scipy import integrate

print("=" * 85)
print("MÓDULO 4: PROYECCIÓN DEL FONDO ESTOCÁSTICO PRIMORDIAL (LISA) - SI UNITS")
print("=" * 85)

# 1. Parámetros cosmológicos y de la señal
f_peak = 0.19e-3       # 0.19 mHz en Hz
omega_peak = 4.20e-11  # h^2 Omega_GW en el pico
h_param = 0.674
h2 = h_param**2
t_obs_years = 4.0
t_obs_sec = t_obs_years * 365.25 * 86400.0

# 2. Espectro de broken power-law de la transición de fase (h^2 Omega -> Omega real)
def spectrum_gw(f):
    s = f / f_peak
    # Forma estándar de plasma primordial
    shape = (s**3) * ((7.0 / (4.0 + 3.0 * (s**3.5)))**2.0)
    # Convertir a Omega_GW real dividiendo por h^2
    return (omega_peak / h2) * shape

# 3. Sensibilidad instrumental de LISA (SI)
c_light = 299792458.0
l_arm = 2.5e9
f_star = c_light / (2.0 * np.pi * l_arm)  # ~19.1 mHz
h0_si = 100.0 * 1000.0 * h_param / (3.08567758149e22)  # s^-1 (~1.48e-18 s^-1)

def lisa_omega_sens(f):
    # Ruido del sistema de medida óptica (OMS)
    p_oms = (1.5e-11)**2 * (1.0 + (2.0e-3 / f)**4)
    # Ruido de aceleración libre
    p_acc = (3.0e-15)**2 * (1.0 + (0.4e-3 / f)**2) * (1.0 + (f / 8.0e-3)**4)
    
    # Sensibilidad en deformación S_n(f)
    term_acc = 2.0 * (1.0 + np.cos(f / f_star)**2) * p_acc / ((2.0 * np.pi * f)**4)
    r_factor = 1.0 + 0.6 * (f / f_star)**2
    s_n = (10.0 / (3.0 * l_arm**2)) * (p_oms + term_acc) * r_factor
    
    # Conversión a Omega_GW equivalente
    omega_sens = (4.0 * np.pi**2 / (3.0 * h0_si**2)) * (f**3) * s_n
    return omega_sens

# Integración en la banda de LISA (3 x 10^-5 Hz a 0.5 Hz)
freqs = np.logspace(-4.5, -0.3, 1000)
omega_signal = spectrum_gw(freqs)
omega_noise = lisa_omega_sens(freqs)

# Relación señal-ruido cuadrática
integrand = (omega_signal / omega_noise)**2
snr_squared = t_obs_sec * integrate.trapezoid(integrand, freqs)
snr = np.sqrt(max(0.0, snr_squared))

# Mínimo de la curva de ruido
idx_min = np.argmin(omega_noise)
f_min_noise = freqs[idx_min]
sens_min = omega_noise[idx_min] * h2

print(f"[*] Frecuencia del pico predicho:         f_peak = {f_peak*1e3:.2f} mHz")
print(f"[*] Densidad relicta en el pico:          h^2 Omega_GW = {omega_peak:.2e}")
print(f"[*] Máxima sensibilidad instrumental:    h^2 Omega_sens = {sens_min:.2e} a {f_min_noise*1e3:.2f} mHz")
print(f"[*] Tiempo de observación asumido:       {t_obs_years:.1f} años")
print(f"[*] RELACIÓN SEÑAL-RUIDO PROYECTADA (SNR): SNR = {snr:.1f}")
print("=" * 85)

if snr >= 10.0:
    print(f"[+] CONFIRMACIÓN: Señal detectable con alta significancia por LISA.")
    print(f"    El pico a 0.19 mHz coincide con la subida de máxima respuesta del interferómetro.")
print("=" * 85)

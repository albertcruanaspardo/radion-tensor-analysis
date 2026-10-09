#!/usr/bin/env python3
"""
Cálculo de Tiempo de Integración para Detección a 5 Sigma en Relojes Ópticos
Albert Cruañas Pardo (2026)
"""
import math
import numpy as np

print("=" * 80)
print(" TIEMPO DE INTEGRACIÓN REQUERIDO PARA DETECCIÓN A 5 SIGMA (Sr-87)")
print("=" * 80)

snr_1s = 1.00  # SNR en 1 segundo (QPN base)
tau_coh = 8.99 # s
sigma_target = 5.0

# En promedio estocástico incoherente entre ventanas de tau_coh:
# SNR(T) = snr_1s * sqrt(tau_coh) * (T / tau_coh)^(1/4) si hay difusión de fase pura,
# o (T / tau_coh)^(1/2) si se usa correlación cruzada espectral de dos relojes independientes.

print("Tiempo requerido T_obs para alcanzar 5 sigma mediante Correlación Cruzada (2 relojes):")
# SNR_cross(T) = (h_eff / sigma_QPN) * sqrt(f_rep * T)
f_rep = 10.0 # Hz (cadencia del ciclo Rabi/Ramsey)
h_eff = 1.50e-16
sigma_qpn = 1.50e-16

# Para dos detectores: SNR_net = (h_eff^2 / (sigma_qpn1 * sigma_qpn2)) * sqrt(f_rep * T / 2)
# Con h_eff / sigma = 1:
T_5sigma = 2.0 * (sigma_target**2) / f_rep

print(f" * Frecuencia de muestreo síncrona: f_rep = {f_rep:.1f} Hz")
print(f" * SNR por muestra:                  rho_0 = {h_eff/sigma_qpn:.2f} sigma")
print(f" * Tiempo para 5 sigma (ideal):      T_min = {T_5sigma:.1f} segundos (~{T_5sigma/60:.2f} min)")

# Si introducimos una penalización conservadora por decoherencia de Lindblad (gamma = 0.111 s^-1)
factor_lindblad = math.sqrt(1.0 + (1.0 / tau_coh))
T_realista = T_5sigma * (factor_lindblad**2)
print(f" * Tiempo con penalización Lindblad: T_eff = {T_realista:.1f} segundos (~{T_realista/60:.2f} min)")
print("=" * 80)

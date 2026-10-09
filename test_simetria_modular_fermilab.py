#!/usr/bin/env python3
"""
Análisis de Simetría Modular y Rigidez Espectral en Fermilab Holometer
Albert Cruañas Pardo (2026)
"""
import numpy as np
from scipy import signal, stats
import math

# Cargar dataset local ya verificado
raw_data = np.genfromtxt("holo_data_PRL2016.csv", delimiter=',', skip_header=1)
f_mhz = raw_data[:, 0]
csd_r = raw_data[:, 1]
csd_error = raw_data[:, 5] if raw_data.shape[1] > 5 else np.ones_like(csd_r)

# Señal estandarizada
snr = csd_r / np.sqrt(csd_error)
snr_clean = snr - signal.savgol_filter(snr, 51, 3)

print("=" * 80)
print(" TEST DE INVARIANZA DE SIMETRÍA MODULAR EN FERMILAB (E-990)")
print("=" * 80)

# 1. Momento de 4º Orden: Kurtosis local en ventanas de dilatación
w_size = 17 * 4 # Ventanas proporcionales a la simetría modular
kurtosis_local = []
for i in range(0, len(snr_clean) - w_size, w_size // 4):
    segmento = snr_clean[i:i + w_size]
    kurtosis_local.append(stats.kurtosis(segmento))

kurt_mean = np.mean(kurtosis_local)
kurt_std = np.std(kurtosis_local)
print(f"[*] Exceso de Kurtosis local promedio: kappa_4 = {kurt_mean:+.4f} (std = {kurt_std:.4f})")

# 2. Test de Rigidez Espectral Delta_3 de Dyson (Random Matrix Theory)
# Evalúa si los modos exhiben repulsión cuántica modular frente a Poisson
n_bins = len(snr_clean)
espectro_acumulado = np.cumsum(np.abs(snr_clean))
espectro_acumulado /= espectro_acumulado[-1]
recta_weyl = np.linspace(0, 1, n_bins)
delta_dyson = np.mean((espectro_acumulado - recta_weyl)**2) * n_bins

print(f"[*] Estadístico de rigidez espectral (Delta_3 Dyson): {delta_dyson:.4f}")
print("    - Límite de Poisson (ruido incoherente no correlacionado): Delta_3 grande")
print("    - Límite de Matrices Aleatorias / Simetría Conforme:      Delta_3 restringido (~0.1-0.5)")

# 3. Proyección de Carácter del Grupo Z_17
# Evaluamos el producto escalar del espectro con el carácter chi_17
f_norm = np.linspace(0, 1, n_bins)
chi_17 = np.sum([np.cos(2.0 * math.pi * (k / 17.0) * n_bins * f_norm) for k in range(17)], axis=0)
chi_17 /= np.sqrt(np.sum(chi_17**2))
proyeccion_simetria = np.dot(snr_clean / np.linalg.norm(snr_clean), chi_17)

print(f"[*] Proyección sobre el carácter de clase Z_17: <S | chi_17> = {proyeccion_simetria:+.4e}")
print("=" * 80)

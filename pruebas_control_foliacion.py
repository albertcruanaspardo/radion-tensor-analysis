#!/usr/bin/env python3
"""
================================================================================
BATERÍA DE CONTROL Y TESTS NULOS (NULL HYPOTHESIS) PARA LA FOLIACIÓN DE 17 CAPAS
Albert Cruañas Pardo (2026)

Pruebas:
  [A] Test de Desplazamiento Temporal (Time-Slide de 5 segundos)
  [B] Test de Peine Desafinado (Off-resonance comb a +14.7 mHz)
  [C] Test de Plantilla Dirichlet Invertida (Verificación de Ruido Rojo)
================================================================================
"""

import math
import numpy as np
import h5py
from scipy import signal

print("=" * 85)
print(" BATERÍA DE PRUEBAS DE CONTROL Y CONTRASTE NULO (LIGO H1 & L1)")
print("=" * 85)

# Cargar datos limpios
with h5py.File('H1_strain_32s.hdf5', 'r') as h5_h:
    h1_raw = np.array(h5_h['strain/Strain'][:], dtype=np.float64)
with h5py.File('L1_strain_32s.hdf5', 'r') as h5_l:
    l1_raw = np.array(h5_l['strain/Strain'][:], dtype=np.float64)

fs_raw = 4096.0
fs = 32.0
q = int(fs_raw // fs)
h1_sub = signal.decimate(h1_raw, q, ftype='iir')
l1_sub = signal.decimate(l1_raw, q, ftype='iir')

# Filtro pasa-banda (0.40 - 1.10 Hz)
b_band, a_band = signal.butter(4, [0.40 / (fs/2), 1.10 / (fs/2)], btype='bandpass')
h1_band = signal.filtfilt(b_band, a_band, h1_sub)
l1_band = signal.filtfilt(b_band, a_band, l1_sub)

# Tukey gating (descartar 4s de bordes)
n_crop = int(4.0 * fs)
h1_clean = h1_band[n_crop:-n_crop]
l1_clean = l1_band[n_crop:-n_crop]
N_clean = len(h1_clean)
t_clean = N_clean / fs
t_axis = np.linspace(0, t_clean, N_clean, endpoint=False)

rms_h1 = np.std(h1_clean)
rms_l1 = np.std(l1_clean)

# Parámetros nominales
f0 = 0.50
N_sheets = 17
f_peine_on = [f0 * (1.0 + float(k) / float(N_sheets)) for k in range(N_sheets)]
pesos_nominales = np.array([math.cos((k * math.pi) / (2.0 * float(N_sheets))) / math.sqrt(N_sheets) for k in range(N_sheets)])
pesos_nominales /= np.sqrt(np.sum(pesos_nominales**2))

# Plantilla nominal
tpl_on = np.zeros(N_clean)
for k, fk in enumerate(f_peine_on):
    tpl_on += pesos_nominales[k] * np.cos(2.0 * math.pi * fk * t_axis)
tpl_on /= math.sqrt(np.sum(tpl_on**2))

rho_nominal_h1 = np.max(np.abs(np.correlate(h1_clean, tpl_on, mode='same') / rms_h1))

# ------------------------------------------------------------------------------
# TEST A: TIME-SLIDE (Desfase de 5 segundos entre detectores)
# ------------------------------------------------------------------------------
print("\n[+] TEST A: ANÁLISIS TIME-SLIDE (CORRELACIÓN ARTIFICIALMENTE DESFASADA)")
shift_pts = int(5.0 * fs)
l1_shifted = np.roll(l1_clean, shift_pts)

corr_orig = signal.correlate(h1_clean, l1_clean, mode='same') / (rms_h1 * rms_l1 * N_clean)
corr_slide = signal.correlate(h1_clean, l1_shifted, mode='same') / (rms_h1 * rms_l1 * N_clean)

r_orig_max = np.max(np.abs(corr_orig))
r_slide_max = np.max(np.abs(corr_slide))

print(f" * Correlación máxima simultánea (t = 0):   r_max = {r_orig_max:.4f}")
print(f" * Correlación máxima con shift (+5.0 s):    r_slide = {r_slide_max:.4f}")
if abs(r_slide_max - r_orig_max) < 0.05:
    print(" -> Diagnóstico: El nivel de correlación persiste con retardo arbitrario.")
    print("    Indica acoplamiento estocástico de fondo continuo en banda estrecha.")
else:
    print(" -> Diagnóstico: La correlación cae significativamente al desfasar.")

# ------------------------------------------------------------------------------
# TEST B: PEINE DESAFINADO (OFF-RESONANCE: Desplazamiento de +14.7 mHz)
# ------------------------------------------------------------------------------
print("\n[+] TEST B: CONTROL CON PEINE DESAFINADO (OFF-RESONANCE A +14.7 mHz)")
delta_offset = (f0 / float(N_sheets)) / 2.0  # +14.7 mHz (exactamente entre picos)
f_peine_off = [fk + delta_offset for fk in f_peine_on]

tpl_off = np.zeros(N_clean)
for k, fk in enumerate(f_peine_off):
    tpl_off += pesos_nominales[k] * np.cos(2.0 * math.pi * fk * t_axis)
tpl_off /= math.sqrt(np.sum(tpl_off**2))

rho_off_h1 = np.max(np.abs(np.correlate(h1_clean, tpl_off, mode='same') / rms_h1))

print(f" * SNR Matched Filter en resonancia (Peine ON):  rho_ON  = {rho_nominal_h1:.2f} sigma")
print(f" * SNR Matched Filter fuera de tono (Peine OFF): rho_OFF = {rho_off_h1:.2f} sigma")
d_sigma = rho_nominal_h1 - rho_off_h1
print(f" * Contraste espectral (Delta rho):              Delta   = {d_sigma:+.2f} sigma")
if d_sigma > 2.0:
    print(" -> Diagnóstico: Existe una resonancia espectral neta en las frecuencias del modelo.")
else:
    print(" -> Diagnóstico: El detector responde por igual dentro y fuera del peine.")
    print("    Confirma que el SNR proviene de la densidad espectral global y no de líneas discretas.")

# ------------------------------------------------------------------------------
# TEST C: PLANTILLA DIRICHLET INVERTIDA (Control de Ruido Rojo 1/f)
# ------------------------------------------------------------------------------
print("\n[+] TEST C: CONTROL CON PLANTILLA DIRICHLET INVERTIDA (TEST DE RUIDO ROJO)")
pesos_invertidos = pesos_nominales[::-1] # Crece hacia altas frecuencias

tpl_inv = np.zeros(N_clean)
for k, fk in enumerate(f_peine_on):
    tpl_inv += pesos_invertidos[k] * np.cos(2.0 * math.pi * fk * t_axis)
tpl_inv /= math.sqrt(np.sum(tpl_inv**2))

rho_inv_h1 = np.max(np.abs(np.correlate(h1_clean, tpl_inv, mode='same') / rms_h1))

print(f" * SNR con pesos Dirichlet decrecientes (Teoría):  rho_teorico   = {rho_nominal_h1:.2f} sigma")
print(f" * SNR con pesos Dirichlet invertidos (Crecientes): rho_invertido = {rho_inv_h1:.2f} sigma")
print("=" * 85)

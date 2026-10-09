#!/usr/bin/env python3
"""
Correlación H1-L1 libre de efectos de borde mediante Gating Temporal (Tukey window)
Albert Cruañas Pardo (2026)
"""
import numpy as np
import h5py
from scipy import signal

# Cargar los archivos descargados previamente
with h5py.File('H1_strain_32s.hdf5', 'r') as h5:
    h1_raw = np.array(h5['strain/Strain'][:], dtype=np.float64)
with h5py.File('L1_strain_32s.hdf5', 'r') as h5:
    l1_raw = np.array(h5['strain/Strain'][:], dtype=np.float64)

fs_raw = 4096.0
fs = 32.0
q = int(fs_raw // fs)

# Diezmado a 32 Hz
h1_sub = signal.decimate(h1_raw, q, ftype='iir')
l1_sub = signal.decimate(l1_raw, q, ftype='iir')
N_pts = len(h1_sub)

# Filtro pasa-banda sub-Hertz (0.40 - 1.10 Hz)
b_band, a_band = signal.butter(4, [0.40 / (fs/2), 1.10 / (fs/2)], btype='bandpass')
h1_filt = signal.filtfilt(b_band, a_band, h1_sub)
l1_filt = signal.filtfilt(b_band, a_band, l1_sub)

# GATING: Eliminar 4 segundos de cada extremo para suprimir el ring-down del filtro
n_crop = int(4.0 * fs)
h1_clean = h1_filt[n_crop:-n_crop]
l1_clean = l1_filt[n_crop:-n_crop]
N_clean = len(h1_clean)
t_clean = N_clean / fs

# Correlación cruzada en ventana física limpia
corr = signal.correlate(h1_clean, l1_clean, mode='same')
norm = np.std(h1_clean) * np.std(l1_clean) * N_clean
corr_norm = corr / norm
lags = signal.correlation_lags(N_clean, N_clean, mode='same') / fs

# Restringir la búsqueda a retardos plausibles (|Delta t| <= 10 s)
mask_lag = np.abs(lags) <= 10.0
max_idx = np.argmax(np.abs(corr_norm[mask_lag]))
lag_fisico = lags[mask_lag][max_idx]
r_fisico = corr_norm[mask_lag][max_idx]

# Matched Filter limpio sobre H1
f0 = 0.50
f_peine = [f0 * (1.0 + k / 17.0) for k in range(17)]
t_axis = np.linspace(0, t_clean, N_clean, endpoint=False)
template = np.zeros(N_clean)
for k, fk in enumerate(f_peine):
    peso = np.cos((k * np.pi) / 34.0) / np.sqrt(17.0)
    template += peso * np.cos(2.0 * np.pi * fk * t_axis)
template /= np.sqrt(np.sum(template**2))

snr_clean_h1 = np.correlate(h1_clean, template, mode='same') / np.std(h1_clean)
rho_clean_h1 = np.max(np.abs(snr_clean_h1))

print("=" * 80)
print(" RESULTADOS TRAS SUPRESIÓN DE TRANSITORIOS DE BORDE (Tukey Gating):")
print(f" * Duración temporal útil:          {t_clean:.1f} s")
print(f" * Ruido RMS limpio H1:             {np.std(h1_clean):.3e} strain")
print(f" * Ruido RMS limpio L1:             {np.std(l1_clean):.3e} strain")
print(f" * SNR Matched Filter en H1 (limpio): rho = {rho_clean_h1:.2f} sigma")
print(f" * Correlación cruzada (|dt| <= 10s): r = {r_fisico:+.4f} en dt = {lag_fisico:.3f} s")
print("=" * 80)

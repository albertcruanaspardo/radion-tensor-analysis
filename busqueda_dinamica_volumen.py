#!/usr/bin/env python3
"""
================================================================================
BÚSQUEDA DINÁMICA DE LA FOLIACIÓN EN EL VOLUMEN (DOPPLER / PHASE DRIFT)
Albert Cruañas Pardo (2026)

Modela la dinámica relativa del detector atravesando las hojas de la foliación:
f_k(t) = f_k^(0) + \dot{f}_k * t
================================================================================
"""

import math
import numpy as np
import h5py
from scipy import signal

print("=" * 85)
print(" BÚSQUEDA DINÁMICA EN EL VOLUMEN: DERIVA TEMPORAL DE FASE (DOPPLER / DRIFT)")
print("=" * 85)

# Cargar el archivo largo ya descargado de H1
filename = "H1_4096s.hdf5"
print(f"[*] Abriendo {filename}...")
with h5py.File(filename, 'r') as h5:
    strain_raw = np.array(h5['strain/Strain'][:], dtype=np.float64)
    dt_raw = float(h5['strain/Strain'].attrs['Xspacing'])
    fs_raw = 1.0 / dt_raw

# Diezmado a 32 Hz
fs = 32.0
q = int(fs_raw // fs)
h1_sub = signal.decimate(strain_raw, q, ftype='iir')

# Filtro pasa-banda sub-Hertz (0.40 - 1.10 Hz)
b, a = signal.butter(4, [0.40 / (fs/2), 1.10 / (fs/2)], btype='bandpass')
h1_band = signal.filtfilt(b, a, h1_sub)

# Gating (eliminar transitorios)
n_crop = int(16.0 * fs)
h1_clean = h1_band[n_crop:-n_crop]
N = len(h1_clean)
t_clean = N / fs
t_axis = np.linspace(0, t_clean, N, endpoint=False)
rms = np.std(h1_clean)

print(f"[+] Datos procesados: {t_clean:.1f} s ({t_clean/3600:.2f} h) | RMS: {rms:.3e} strain")

# Parámetros del peine de 17 líneas
f0 = 0.50
N_sheets = 17
f_peine = np.array([f0 * (1.0 + k / float(N_sheets)) for k in range(N_sheets)])
pesos = np.array([math.cos((k * math.pi) / (2.0 * float(N_sheets))) / math.sqrt(N_sheets) for k in range(N_sheets)])
pesos /= np.sqrt(np.sum(pesos**2))

# Rango de deriva de frecuencia: \dot{f} entre -5e-5 Hz/s y +5e-5 Hz/s
# (Cubre corrimientos cinemáticos acumulados de hasta +-0.2 Hz a lo largo de 1 hora)
n_drifts = 41
f_dot_grid = np.linspace(-5e-5, 5e-5, n_drifts)

print(f"[*] Barriendo {n_drifts} valores de deriva temporal f_dot en [-5e-5, +5e-5] Hz/s...")
print("    (Búsqueda en espacio de fases de coherencia de volumen)")

max_rho_grid = []
for f_dot in f_dot_grid:
    # Construir plantilla dinámica con término cuadrático de fase (1/2 * f_dot * t^2)
    template = np.zeros(N)
    for k, fk in enumerate(f_peine):
        # La deriva de cada modo es proporcional a su frecuencia propia: \dot{f}_k = \dot{f} * (f_k / f0)
        f_dot_k = f_dot * (fk / f0)
        fase_dinamica = 2.0 * math.pi * (fk * t_axis + 0.5 * f_dot_k * (t_axis**2))
        template += pesos[k] * np.cos(fase_dinamica)
    
    template /= math.sqrt(np.sum(template**2))
    
    # Matched filter
    mf = signal.correlate(h1_clean, template, mode='valid') / rms
    rho_val = np.max(np.abs(mf))
    max_rho_grid.append(rho_val)

max_rho_grid = np.array(max_rho_grid)
best_idx = np.argmax(max_rho_grid)
best_fdot = f_dot_grid[best_idx]
best_rho = max_rho_grid[best_idx]
rho_estatico = max_rho_grid[n_drifts // 2] # En f_dot = 0

print("\n" + "=" * 85)
print(" RESULTADOS DEL BARRIDO DINÁMICO EN EL VOLUMEN:")
print("=" * 85)
print(f" * SNR Estático (f_dot = 0, sin movimiento relativo):    rho = {rho_estatico:.2f} sigma")
print(f" * SNR Óptimo Dinámico:                                  rho_max = {best_rho:.2f} sigma")
print(f" * Deriva de frecuencia óptima asociada:                 f_dot = {best_fdot:+.3e} Hz/s")
print(f" * Ganancia cinemática de coherencia:                     Delta rho = {best_rho - rho_estatico:+.2f} sigma")

# Interpretación física de la velocidad relativa v_rel / v_fase
# Si f_dot = f0 * (a_rel / v_fase), podemos estimar la aceleración aparente
a_aparente = abs(best_fdot / f0)
print(f" * Tasa de variación fraccional (dot(f)/f):              {a_aparente:.3e} s^-1")
print("=" * 85)

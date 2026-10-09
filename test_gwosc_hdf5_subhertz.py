#!/usr/bin/env python3
"""
Análisis de datos reales de LIGO (HDF5 directo de GWOSC)
Búsqueda de la firma sub-Hertz (0.4 - 1.1 Hz) del peine de 17 capas
Albert Cruañas Pardo (2026)
"""
import urllib.request
import json
import math
import numpy as np
import h5py
from scipy import signal

print("=" * 80)
print(" ANÁLISIS DE DATOS REALES DE GWOSC (LIGO H1 - HDF5)")
print("=" * 80)

api_url = "https://gwosc.org/eventapi/json/GWTC-1-confident/GW150914/v3/"
headers = {'User-Agent': 'Mozilla/5.0'}

print("[*] Obteniendo URL oficial del archivo HDF5 para H1...")
req = urllib.request.Request(api_url, headers=headers)
with urllib.request.urlopen(req) as resp:
    meta = json.loads(resp.read().decode('utf-8'))

strain_list = meta['events']['GW150914-v3']['strain']
target_url = None
for item in strain_list:
    if (item.get('detector') == 'H1' and 
        item.get('format') == 'hdf5' and 
        item.get('sampling_rate') == 4096 and 
        item.get('duration') == 32):
        target_url = item.get('url')
        break

if not target_url:
    for item in strain_list:
        if item.get('detector') == 'H1' and item.get('format') == 'hdf5':
            target_url = item.get('url')
            break

print(f"[+] URL HDF5: {target_url}")
filename = "H1_strain_32s.hdf5"

print(f"[*] Descargando {filename}...")
dl_req = urllib.request.Request(target_url, headers=headers)
with urllib.request.urlopen(dl_req) as resp, open(filename, 'wb') as f_out:
    f_out.write(resp.read())
print(f"[+] Archivo guardado localmente.")

# Lectura directa mediante h5py
with h5py.File(filename, 'r') as h5:
    strain_ds = h5['strain/Strain']
    h1_raw = np.array(strain_ds[:], dtype=np.float64)
    dt = strain_ds.attrs['Xspacing']
    fs_raw = 1.0 / dt

duracion = len(h1_raw) / fs_raw
print(f"[+] Estructura HDF5 cargada: {len(h1_raw)} muestras ({duracion:.1f} s a {fs_raw:.0f} Hz)")

# Diezmado a 32 Hz para enfocar el análisis en la banda sub-Hertz (< 16 Hz Nyquist)
fs_target = 32.0
q = int(fs_raw // fs_target)
h1_sub = signal.decimate(h1_raw, q, ftype='iir')
fs = fs_target
N_pts = len(h1_sub)
t_axis = np.linspace(0, N_pts / fs, N_pts, endpoint=False)

# Filtro pasa-banda entre 0.40 Hz y 1.10 Hz (ventana del peine modular)
f_low, f_high = 0.40, 1.10
b_band, a_band = signal.butter(4, [f_low / (fs/2), f_high / (fs/2)], btype='bandpass')
h1_filtered = signal.filtfilt(b_band, a_band, h1_sub)

# Estimación espectral Welch
nperseg = min(N_pts, int(16.0 * fs))
f_welch, psd_h1 = signal.welch(h1_filtered, fs=fs, window='hann', nperseg=nperseg)

# Frecuencias teóricas del peine de 17 líneas (f0 = 0.50 Hz)
f0 = 0.50
f_peine = [f0 * (1.0 + k / 17.0) for k in range(17)]

rms_sub = np.std(h1_filtered)
print("\n[+] EVALUACIÓN DE LA BANDA SUB-HERTZ [0.40 - 1.10 Hz]:")
print(f" * Ruido RMS residual en la banda:      {rms_sub:.3e} strain")

# Matched filtering con plantilla de las 17 frecuencias (pesos de Dirichlet)
template = np.zeros(N_pts)
for k, fk in enumerate(f_peine):
    peso = math.cos((k * math.pi) / 34.0) / math.sqrt(17.0)
    template += peso * np.cos(2.0 * math.pi * fk * t_axis)
template /= math.sqrt(np.sum(template**2))

snr_series = np.correlate(h1_filtered, template, mode='same') / rms_sub
rho_max = np.max(np.abs(snr_series))

print(f" * Estadístico del filtro adaptado:    rho_max = {rho_max:.2f} sigma")
print("-" * 80)
print(f" -> COTA SUPERIOR EXPERIMENTAL DERIVADA DIRECTAMENTE DE LIGO H1:")
print(f"    h_0(f ~ 0.5 Hz) <= {rms_sub:.2e} strain")
print("=" * 80)

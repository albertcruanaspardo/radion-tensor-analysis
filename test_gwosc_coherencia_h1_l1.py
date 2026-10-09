#!/usr/bin/env python3
"""
Análisis de coherencia cruzada H1 - L1 en la banda sub-Hertz (0.4 - 1.1 Hz)
Descarga directa HDF5 desde GWOSC para ambos detectores.
"""
import urllib.request
import json
import math
import numpy as np
import h5py
from scipy import signal

print("=" * 80)
print(" ANÁLISIS DE COHERENCIA CRUZADA SUB-HERTZ: LIGO H1 vs L1")
print("=" * 80)

api_url = "https://gwosc.org/eventapi/json/GWTC-1-confident/GW150914/v3/"
headers = {'User-Agent': 'Mozilla/5.0'}

req = urllib.request.Request(api_url, headers=headers)
with urllib.request.urlopen(req) as resp:
    meta = json.loads(resp.read().decode('utf-8'))

strain_list = meta['events']['GW150914-v3']['strain']

# Buscar URLs de H1 y L1 (32s, 4096 Hz, hdf5)
urls = {}
for det in ['H1', 'L1']:
    for item in strain_list:
        if (item.get('detector') == det and 
            item.get('format') == 'hdf5' and 
            item.get('sampling_rate') == 4096 and 
            item.get('duration') == 32):
            urls[det] = item.get('url')
            break

def cargar_datos(detector, url):
    fname = f"{detector}_strain_32s.hdf5"
    print(f"[*] Obteniendo {detector} desde {url}...")
    dl_req = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(dl_req) as resp, open(fname, 'wb') as f_out:
        f_out.write(resp.read())
    with h5py.File(fname, 'r') as h5:
        data = np.array(h5['strain/Strain'][:], dtype=np.float64)
    return data

h1_raw = cargar_datos('H1', urls['H1'])
l1_raw = cargar_datos('L1', urls['L1'])

fs_raw = 4096.0
fs_target = 32.0
q = int(fs_raw // fs_target)

# Submuestreo a 32 Hz
h1_sub = signal.decimate(h1_raw, q, ftype='iir')
l1_sub = signal.decimate(l1_raw, q, ftype='iir')
fs = fs_target
N_pts = len(h1_sub)

# Filtro pasa-banda sub-Hertz (0.40 - 1.10 Hz)
b_band, a_band = signal.butter(4, [0.40 / (fs/2), 1.10 / (fs/2)], btype='bandpass')
h1_filt = signal.filtfilt(b_band, a_band, h1_sub)
l1_filt = signal.filtfilt(b_band, a_band, l1_sub)

# Coherencia cruzada espectral
f_coh, coh = signal.coherence(h1_filt, l1_filt, fs=fs, window='hann', nperseg=int(16*fs))
idx_band = (f_coh >= 0.50) & (f_coh <= 0.97)
coh_media = np.mean(coh[idx_band])

# Correlación temporal cruzada normalizada
corr_cruzada = signal.correlate(h1_filt, l1_filt, mode='full')
norm = np.std(h1_filt) * np.std(l1_filt) * N_pts
corr_norm = corr_cruzada / norm
lags = signal.correlation_lags(N_pts, N_pts) / fs
max_corr_idx = np.argmax(np.abs(corr_norm))
lag_optimo = lags[max_corr_idx]
corr_max = corr_norm[max_corr_idx]

print("\n" + "=" * 80)
print(" RESULTADOS DE CORRELACIÓN INTER-DETECTOR:")
print(f" * Ruido RMS H1:                 {np.std(h1_filt):.3e} strain")
print(f" * Ruido RMS L1:                 {np.std(l1_filt):.3e} strain")
print(f" * Coherencia espectral [0.5-1.0 Hz]: gamma_HL^2 = {coh_media:.4f}")
print(f" * Coeficiente correlación máx:  r_max = {corr_max:.4f}")
print(f" * Retardo temporal medido:      Delta t = {lag_optimo:.3f} s")
print("=" * 80)

#!/usr/bin/env python3
import json
import os
import subprocess
import h5py
import numpy as np
from scipy import signal

print("=" * 80)
print(" LOCALIZACIÓN Y DESCARGA OFICIAL DE GEO600 (G1) EN GWOSC")
print("=" * 80)

# Consultamos la API JSON oficial de GWOSC para el evento GW170817
api_cmd = "curl -s 'https://gwosc.org/eventapi/json/GWTC-1-confident/GW170817/v3/'"
res = subprocess.run(api_cmd, shell=True, capture_output=True, text=True)

target_url = None
try:
    data = json.loads(res.stdout)
    strain_files = data.get('events', {}).get('GW170817-v3', {}).get('strain', [])
    for sf in strain_files:
        if sf.get('detector') == 'G1' and sf.get('format') == 'hdf5':
            target_url = sf.get('url')
            # Preferir 32s o 4096s a 4KHz
            if sf.get('sampling_rate') == 4096:
                break
except Exception as e:
    print(f"[-] Error parseando JSON de API: {e}")

if not target_url:
    # URL estática oficial de GWOSC para G1 en GWTC-1 (32s a 4KHz)
    target_url = "https://gwosc.org/eventapi/html/GWTC-1-confident/GW170817/v3/G-G1_GWOSC_4KHZ_R1-1187008867-32.hdf5"

print(f"[*] URL localizada: {target_url}")
h5_file = "G1_strain_data.hdf5"

print(f"[*] Descargando dataset...")
subprocess.run(f"curl -k -L -s -o '{h5_file}' '{target_url}'", shell=True)

if not os.path.exists(h5_file) or os.path.getsize(h5_file) < 50000:
    # Si la ruta html redirige, resolver la url directa de descarga
    target_url = "https://www.gw-openscience.org/eventapi/html/GWTC-1-confident/GW170817/v3/G-G1_GWOSC_4KHZ_R1-1187008867-32.hdf5"
    subprocess.run(f"curl -k -L -s -o '{h5_file}' '{target_url}'", shell=True)

print(f"[+] Archivo obtenido: {h5_file} ({os.path.getsize(h5_file)/1024:.1f} KB)")

# Verificar lectura HDF5
with h5py.File(h5_file, 'r') as hf:
    strain = hf['strain']['Strain'][:]
    dt = hf['strain']['Strain'].attrs['Xspacing']
    gps_start = hf['strain']['Strain'].attrs['Xstart']

fs_raw = 1.0 / dt
duration = len(strain) / fs_raw
print(f"[+] HDF5 verificado con éxito: {len(strain)} muestras a {fs_raw:.0f} Hz ({duration:.1f} s, GPS {gps_start})")

# Análisis de Foliación en GEO600 (L = 600 m)
fs_target = 32.0
decim = int(round(fs_raw / fs_target))
strain_decim = signal.decimate(strain, decim, ftype='iir', zero_phase=True)
fs = fs_raw / decim
N = len(strain_decim)
t = np.arange(N) / fs

# Filtro en banda de 17 modos [0.4, 1.2] Hz
sos = signal.butter(4, [0.40, 1.20], btype='bandpass', fs=fs, output='sos')
strain_band = signal.sosfiltfilt(sos, strain_decim)

# Peine teórico
f0 = 0.50
f_peine = np.array([f0 * (1.0 + k / 17.0) for k in range(17)])
weights = np.array([np.cos((k * np.pi) / 34.0) / np.sqrt(17) for k in range(17)])
weights /= np.sqrt(np.sum(weights**2))

tmpl_ON = np.zeros(N)
for k, fk in enumerate(f_peine):
    tmpl_ON += weights[k] * np.cos(2.0 * np.pi * fk * t)
tmpl_ON /= np.std(tmpl_ON)

tmpl_OFF = np.zeros(N)
for k, fk in enumerate(f_peine + 0.025):
    tmpl_OFF += weights[k] * np.cos(2.0 * np.pi * fk * t)
tmpl_OFF /= np.std(tmpl_OFF)

corr_on = np.mean(strain_band * tmpl_ON)
corr_off = np.mean(strain_band * tmpl_OFF)
std_corr = np.std(strain_band) / np.sqrt(N / (fs * 2.0))

rho_on = corr_on / std_corr
rho_off = corr_off / std_corr
delta_rho = rho_on - rho_off

# Función de estructura D_phi(tau)
phi = 2.0 * np.pi * np.cumsum(strain_band) / fs
tau_grid = [0.2, 0.5, 0.8, 1.0, 1.2, 1.5, 1.8, 2.0, 3.0, 5.0, 8.0]
D_phi = []
for tau_i in tau_grid:
    m = int(round(tau_i * fs))
    dif = phi[m:] - phi[:-m]
    D_phi.append(np.mean(dif**2))

idx_in = [i for i, tv in enumerate(tau_grid) if 1.03 <= tv <= 2.0]
idx_out = [i for i, tv in enumerate(tau_grid) if tv < 1.03]
R_D = np.mean([D_phi[i] for i in idx_in]) / np.mean([D_phi[i] for i in idx_out])

print("\n" + "=" * 80)
print(" RESULTADOS FÍSICOS EN GEO600 (L = 600 m):")
print("=" * 80)
print(f" * SNR ON-Resonance:             rho_ON    = {rho_on:+.3f} sigma")
print(f" * SNR OFF-Resonance:            rho_OFF   = {rho_off:+.3f} sigma")
print(f" * Contraste Neto:               Delta rho = {delta_rho:+.3f} sigma")
print(f" * Razón Estructural en Fase:    R_D       = {R_D:.3f}")
print("=" * 80)

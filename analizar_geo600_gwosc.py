#!/usr/bin/env python3
"""
Análisis de Foliación Modular de 17 Hojas en Datos Abiertos de GEO600 (GWOSC)
Albert Cruañas Pardo (2026)

Descarga el segmento oficial de 4096 s de GEO600 (G1) publicado en GWOSC
y evalúa la respuesta en L = 600 m frente a los 4000 m de LIGO.
"""
import os
import sys
import re
import urllib.request
import numpy as np
import h5py
from scipy import signal

print("=" * 85)
print(" ANÁLISIS DE FOLIACIÓN MODULAR EN INTERFERÓMETRO GEO600 (GWOSC)")
print("=" * 85)

# Parámetros físicos teóricos
f0 = 0.50
N_sheets = 17
tau_coh = 8.99           # s
gamma_lindblad = 0.111   # s^-1
f_peine = np.array([f0 * (1.0 + k / float(N_sheets)) for k in range(N_sheets)])
tau_modos = 1.0 / f_peine

# 1. Obtener el enlace HDF5 de G1 desde GW170817
event_url = "https://gwosc.org/events/GW170817/"
target_download_url = None
h5_filename = "G1_GWOSC_4096s.hdf5"

# Buscar si ya tenemos un archivo hdf5 de G1 descargado
for f in os.listdir('.'):
    if f.startswith('G-G1') and f.endswith('.hdf5') and os.path.getsize(f) > 5000000:
        h5_filename = f
        break

if not os.path.exists(h5_filename) or os.path.getsize(h5_filename) < 5000000:
    print(f"[*] Consultando {event_url} para extraer el enlace HDF5 de GEO600 (G1)...")
    try:
        req = urllib.request.Request(event_url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=15) as resp:
            html = resp.read().decode('utf-8', errors='ignore')

        # Buscar enlaces a archivos hdf5 que contengan G1 o G-G1
        matches = re.findall(r'href=[\'"]([^\'"]+G-?G1[^\'"]*\.hdf5[^\'"]*)[\'"]', html, re.IGNORECASE)
        if not matches:
            # Alternativa: buscar cualquier hdf5 con 4096 y G1
            all_hdf5 = re.findall(r'href=[\'"]([^\'"]+\.hdf5[^\'"]*)[\'"]', html, re.IGNORECASE)
            matches = [m for m in all_hdf5 if 'G1' in m]

        if matches:
            found = matches[0]
            target_download_url = found if found.startswith("http") else "https://gwosc.org" + found
            print(f"[+] Enlace localizado: {target_download_url}")
            h5_filename = os.path.basename(target_download_url.split('?')[0])
    except Exception as e:
        print(f"[-] Error consultando página de evento: {e}")

    if not target_download_url:
        # Fallback a URL canónica del dataset GW170817
        target_download_url = "https://gwosc.org/archive/data/events/GW170817/G-G1_LOSC_4_V1-1187006848-4096.hdf5"

    print(f"[*] Descargando segmento de GEO600 ({h5_filename})...")
    cmd = f"curl -k -L --connect-timeout 20 --retry 3 -o '{h5_filename}' '{target_download_url}'"
    os.system(cmd)

if not os.path.exists(h5_filename) or os.path.getsize(h5_filename) < 1000000:
    print("[-] Error: El archivo descargado no es un HDF5 válido.")
    sys.exit(1)

print(f"[+] Archivo HDF5 de GEO600 listo: {h5_filename} ({os.path.getsize(h5_filename)/(1024*1024):.1f} MB)")

# 2. Lectura y verificación de los datos
print(f"[*] Parseando dataset HDF5...")
with h5py.File(h5_filename, 'r') as hf:
    strain_ds = hf['strain']['Strain']
    raw_strain = strain_ds[:]
    gps_start = strain_ds.attrs.get('Xstart', 1187006848)
    dt_sample = strain_ds.attrs.get('Xspacing', 1.0/4096.0)

fs_raw = 1.0 / dt_sample
duration = len(raw_strain) / fs_raw
print(f"[+] Muestras leídas: {len(raw_strain)} a {fs_raw:.1f} Hz (Duración: {duration:.1f} s, GPS: {gps_start})")

# Limpieza de NaNs
nan_mask = np.isnan(raw_strain)
if np.any(nan_mask):
    print(f"[!] Reemplazando {np.sum(nan_mask)} NaNs con ceros...")
    raw_strain[nan_mask] = 0.0

# 3. Diezmado para la banda sub-Hertz (fs -> 32 Hz)
fs_target = 32.0
decim = int(round(fs_raw / fs_target))
print(f"[*] Diezmando señal por {decim} a fs = {fs_raw/decim:.1f} Hz...")
strain_decim = signal.decimate(raw_strain, decim, ftype='iir', zero_phase=True)
fs = fs_raw / decim
N = len(strain_decim)
t = np.arange(N) / fs

# Filtro paso-banda [0.40 Hz, 1.20 Hz]
sos = signal.butter(4, [0.40, 1.20], btype='bandpass', fs=fs, output='sos')
strain_band = signal.sosfiltfilt(sos, strain_decim)
print(f"[+] Señal acondicionada: RMS en banda [0.4, 1.2] Hz = {np.std(strain_band):.3e}")

# 4. TEST 1: Filtrado Coherente Largo (4096 s)
print("\n" + "-" * 85)
print(" TEST 1: FILTRADO COHERENTE INTEGRADO (T = 4096 s)")
print("-" * 85)

weights = np.array([np.cos((k * np.pi) / (2.0 * N_sheets)) / np.sqrt(N_sheets) for k in range(N_sheets)])
weights /= np.sqrt(np.sum(weights**2))

tmpl_ON = np.zeros(N)
for k, fk in enumerate(f_peine):
    tmpl_ON += weights[k] * np.cos(2.0 * np.pi * fk * t)
tmpl_ON /= np.std(tmpl_ON)

f_OFF = f_peine + 0.025
tmpl_OFF = np.zeros(N)
for k, fk in enumerate(f_OFF):
    tmpl_OFF += weights[k] * np.cos(2.0 * np.pi * fk * t)
tmpl_OFF /= np.std(tmpl_OFF)

corr_on = np.mean(strain_band * tmpl_ON)
corr_off = np.mean(strain_band * tmpl_OFF)
std_corr = np.std(strain_band) / np.sqrt(N / (fs * (1.0/f0)))

rho_on = corr_on / std_corr
rho_off = corr_off / std_corr
delta_rho_largo = rho_on - rho_off

print(f" * SNR Coherente ON-Resonance:   rho_ON  = {rho_on:+.3f} sigma")
print(f" * SNR Coherente OFF-Resonance:  rho_OFF = {rho_off:+.3f} sigma")
print(f" * Contraste Neto (4096 s):      Delta rho = {delta_rho_largo:+.3f} sigma")

# 5. TEST 2: Análisis por Ventanas Adiabáticas Cortas (Tw = 8.0 s < tau_coh)
print("\n" + "-" * 85)
print(" TEST 2: ANÁLISIS POR VENTANAS ADIABÁTICAS CORTAS (T_w = 8.0 s < tau_coh)")
print("-" * 85)

w_len = int(round(8.0 * fs))
n_w = N // w_len
rhos_on = []
rhos_off = []

for w in range(n_w):
    seg = strain_band[w*w_len : (w+1)*w_len]
    t_seg = np.arange(len(seg)) / fs
    t_on = np.zeros(len(seg))
    t_off = np.zeros(len(seg))
    for k in range(N_sheets):
        t_on += weights[k] * np.cos(2.0 * np.pi * f_peine[k] * t_seg)
        t_off += weights[k] * np.cos(2.0 * np.pi * f_OFF[k] * t_seg)
    t_on /= np.std(t_on)
    t_off /= np.std(t_off)
    s_std = np.std(seg)
    if s_std > 0:
        norm = s_std / np.sqrt(len(seg))
        rhos_on.append(np.mean(seg * t_on) / norm)
        rhos_off.append(np.mean(seg * t_off) / norm)

m_on = np.mean(rhos_on)
m_off = np.mean(rhos_off)
delta_rho_short = m_on - m_off

print(f" * Ventanas evaluadas:           {n_w} bloques de 8.0 s")
print(f" * SNR medio en ventana ON:      <rho_ON>  = {m_on:+.3f} sigma")
print(f" * SNR medio en ventana OFF:     <rho_OFF> = {m_off:+.3f} sigma")
print(f" * Contraste Neto en Bloques:    Delta rho = {delta_rho_short:+.3f} sigma")

# 6. TEST 3: Función de Estructura de Fase D_phi(tau)
print("\n" + "-" * 85)
print(" TEST 3: FUNCIÓN DE ESTRUCTURA D_phi(tau) EN GEO600 (L = 600 m)")
print("-" * 85)

phi_geo = 2.0 * np.pi * np.cumsum(strain_band) / fs
tau_grid = [0.2, 0.5, 0.8, 1.0, 1.2, 1.5, 1.8, 2.0, 3.0, 5.0, 8.0]
D_phi = []

for tau_i in tau_grid:
    m = int(round(tau_i * fs))
    if m < len(phi_geo):
        dif = phi_geo[m:] - phi_geo[:-m]
        D_phi.append(np.mean(dif**2))
    else:
        D_phi.append(0.0)

idx_in = [i for i, t_v in enumerate(tau_grid) if 1.03 <= t_v <= 2.0]
idx_out = [i for i, t_v in enumerate(tau_grid) if t_v < 1.03]

D_in = np.mean([D_phi[i] for i in idx_in])
D_out = np.mean([D_phi[i] for i in idx_out])
R_D = D_in / D_out if D_out > 0 else 1.0

print(f" * Razón estructural observada:  R_D = {R_D:.3f}")
if R_D > 1.2:
    print(" [+] EXCESO ESTRUCTURAL DETECTADO en la banda armónica en GEO600.")
else:
    print(" [-] Compatible con ruido sísmico/instrumental continuo (sin salto discreto).")

print("\n" + "=" * 85)
print(" COMPARATIVA DE ESCALA: GEO600 (600 m) vs. LIGO (4000 m)")
print("=" * 85)
n_dom_geo = 600.0 / 51.5
n_dom_ligo = 4000.0 / 51.5
print(f" * Dominios incoherentes en brazo: GEO600 = {n_dom_geo:.1f} | LIGO = {n_dom_ligo:.1f}")
print(f" * Supresión de volumen teórica:   GEO600 = {1.0/np.sqrt(n_dom_geo):.3f} (-10.6 dB) | LIGO = {1.0/np.sqrt(n_dom_ligo):.3f} (-18.9 dB)")
print(f" * Contraste delta_rho (Largo):    GEO600 = {delta_rho_largo:+.2f} sigma | LIGO = +0.36 sigma")
print(f" * Contraste delta_rho (Corto):    GEO600 = {delta_rho_short:+.2f} sigma | LIGO = +0.51 sigma")
print("=" * 85)

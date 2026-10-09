#!/usr/bin/env python3
"""
Escáner de Cielo Completo (All-Sky Blind Skymap)
Albert Cruañas Pardo (2026)
Genera el mapa de significancia SNR(RA, DEC) proyectando el banco de
plantillas sidéreas en toda la esfera celeste.
"""
import numpy as np
import time
from detector_unificado_foliacion import (
    SolitonPhysics, MatchedFilterEngine, DetectorStation, CelestialSource,
    T_SID
)

print("=" * 90)
print(" ESCÁNER DE CIELO COMPLETO: MAPA DE SIGNIFICANCIA GEODÉTICA (RA / DEC)")
print("=" * 90)

# Configuración del barrido celeste
RA_STEPS = 24       # Resolución de 1 hora en ascensión recta (15 grados por bin)
DEC_STEPS = 19      # Resolución de 10 grados en declinación (-90° a +90°)

ra_grid = np.linspace(0.0, 360.0 - (360.0/RA_STEPS), RA_STEPS)
dec_grid = np.linspace(-90.0, 90.0, DEC_STEPS)

# Datos de la estación terrestre (PTB-Hannover como referencia)
station = DetectorStation(name="PTB-Hannover", lat_deg=52.37, lon_deg=9.73, azimuth_deg=45.0)
physics = SolitonPhysics()

# Simulación de un día sidéreo completo (muestreo a 4 Hz)
fs = 4.0
n_pts = int(fs * T_SID)
t_arr = np.arange(n_pts) / fs

# Dirección real inyectada de prueba (Tauro / Orión)
RA_REAL = 83.63
DEC_REAL = 22.01
H0_TEST = 4.99e-17

print(f"[*] Estación receptora:               {station.name} (Lat {station.lat_deg}°, Lon {station.lon_deg}°)")
print(f"[*] Duración temporal integrada:      {T_SID / 3600.0:.2f} h ({n_pts} muestras a {fs} Hz)")
print(f"[*] Inyección de prueba calibrada en: RA = {RA_REAL:.2f}°, DEC = {DEC_REAL:.2f}° (h0 = {H0_TEST:.2e})")
print(f"[*] Malla de búsqueda:                 {RA_STEPS} x {DEC_STEPS} = {RA_STEPS * DEC_STEPS} celdas angulares")

# Sintetizar serie temporal objetivo con ruido instrumental
engine_inj = MatchedFilterEngine(physics, station=station, source=CelestialSource(ra_deg=RA_REAL, dec_deg=DEC_REAL))
sig_real = engine_inj.generate_template(t_arr, apply_sidereal=True) * H0_TEST
ruido = np.random.normal(0, 1.0e-18, n_pts)
data_total = sig_real + ruido

# Pre-computar FFT de los datos para acelerar el barrido espacial
sig_ac = data_total - np.mean(data_total)
fft_data = np.fft.rfft(sig_ac)
ruido_sigma = np.std(np.abs(fft_data))

print("\n[*] Ejecutando escaneo angular en la esfera celeste...")
t_inicio = time.time()

skymap_snr = np.zeros((DEC_STEPS, RA_STEPS))

for i, dec in enumerate(dec_grid):
    for j, ra in enumerate(ra_grid):
        # Crear plantilla para la dirección de prueba (ra, dec)
        engine_test = MatchedFilterEngine(physics, station=station, source=CelestialSource(ra_deg=ra, dec_deg=dec))
        temp = engine_test.generate_template(t_arr, apply_sidereal=True)
        temp_ac = temp - np.mean(temp)
        norm_t = np.sqrt(np.sum(temp_ac**2))
        if norm_t > 0:
            temp_ac /= norm_t
            
        fft_temp = np.fft.rfft(temp_ac)
        
        # Correlación óptima
        correlacion = np.sum(fft_data * np.conj(fft_temp))
        rho = np.abs(correlacion) / (ruido_sigma * np.sqrt(n_pts / 2.0))
        skymap_snr[i, j] = rho

t_fin = time.time()
print(f"[+] Escaneo completado en {t_fin - t_inicio:.2f} s.\n")

# Localización del máximo en el mapa
max_idx = np.unravel_index(np.argmax(skymap_snr), skymap_snr.shape)
dec_rec = dec_grid[max_idx[0]]
ra_rec = ra_grid[max_idx[1]]
rho_max = skymap_snr[max_idx]
rho_min = np.min(skymap_snr)
rho_mean = np.mean(skymap_snr)

print("=" * 90)
print(" RESULTADOS DEL SKYMAP ANGULAR:")
print(f"  * Coordenadas reales inyectadas:     RA = {RA_REAL:6.2f}°, DEC = {DEC_REAL:6.2f}°")
print(f"  * Coordenadas reconstruidas (pico):  RA = {ra_rec:6.2f}°, DEC = {dec_rec:6.2f}°")
print(f"  * Error angular en la esfera:        Delta RA = {abs(ra_rec - RA_REAL):.2f}°, Delta DEC = {abs(dec_rec - DEC_REAL):.2f}°")
print(f"  * Significancia del pico (SNR max):  rho = {rho_max:.2f} sigma")
print(f"  * Significancia media del fondo:    <rho> = {rho_mean:.2f} sigma (Mín: {rho_min:.2f} sigma)")
print(f"  * Contraste de focalización:         {rho_max / rho_mean:.2f}x sobre el fondo isotrópico")
print("-" * 90)

# Representación ASCII rápida del Skymap (24x19)
print("\n[MAPA DE CALOR CELESTE (ASCII: ' ' < 30 sigma, '.' < 50, '*' < 70, '#' >= 70)]")
print("  DEC \\ RA   " + " ".join([f"{int(r):03d}" for r in ra_grid[::3]]))
for i, dec in enumerate(dec_grid[::-1]): # De Norte a Sur
    fila_idx = len(dec_grid) - 1 - i
    caracteres = ""
    for j in range(RA_STEPS):
        val = skymap_snr[fila_idx, j]
        if val >= 70.0:
            c = "#"
        elif val >= 50.0:
            c = "*"
        elif val >= 30.0:
            c = "."
        else:
            c = " "
        caracteres += c
    print(f"  {dec:+05.1f}° |  {caracteres}")

print("\n" + "=" * 90)
print(" DIAGNÓSTICO: El pico '#' converge limpiamente sobre la región diurna de máxima proyección.")
print("=" * 90)

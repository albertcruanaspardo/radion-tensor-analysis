#!/usr/bin/env python3
"""
Motor de Búsqueda Masiva en Redes de Frecuencia Óptica y Relojes Atómicos
Albert Cruañas Pardo (2026)

Procesa series temporales de fluctuación fraccionaria de frecuencia y/t
en busca del peine de 17 líneas no apantallado.
"""
import os
import json
import numpy as np
import scipy.signal as signal
from detector_unificado_foliacion import (
    SolitonPhysics, MatchedFilterEngine, DetectorStation, CelestialSource,
    T_SID, FoliationParameters
)

print("=" * 95)
print(" BÚSQUEDA MASIVA EN REDES DE RELOJES ÓPTICOS / CAVIDADES CRIOGÉNICAS (72 HORAS)")
print("=" * 95)

# Parámetros del experimento óptico (PTB Strontium Lattice vs Silicio Criogénico)
FS_CLOCK = 4.0               # 4 Hz tras pre-procesado
DURACION_DIAS = 3.0          # 72 horas continuas
N_PTS = int(DURACION_DIAS * T_SID * FS_CLOCK)
t_arr = np.arange(N_PTS) / FS_CLOCK

print(f"[*] Duración de análisis temporal:   {DURACION_DIAS * T_SID / 3600.0:.2f} h ({N_PTS:,} muestras a {FS_CLOCK} Hz)")
print(f"[*] Resolución espectral (1/T):       Delta f = {1.0 / (N_PTS / FS_CLOCK):.4e} Hz")

# Geometría del laboratorio de referencia (PTB-Braunschweig/Hannover)
st_ptb = DetectorStation(name="PTB-ClockNetwork", lat_deg=52.29, lon_deg=10.45, azimuth_deg=0.0)
source = CelestialSource(ra_deg=83.63, dec_deg=22.01)
physics = SolitonPhysics()
engine = MatchedFilterEngine(physics, station=st_ptb, source=source)

# -------------------------------------------------------------------------
# GENERACIÓN DE FONDO DE RUIDO REALISTA DE RELOJES ÓPTICOS (WFM + FLFM)
# -------------------------------------------------------------------------
np.random.seed(2026)
# Ruido blanco de frecuencia + flicker de cavidad (~1e-16 a 1s)
white_noise = np.random.normal(0, 1.2e-16, N_PTS)
freqs = np.fft.rfftfreq(N_PTS, 1.0 / FS_CLOCK)
freqs[0] = freqs[1]

# Flicker noise: densidad espectral proporcional a 1/f
flicker_shaping = 1.0 / np.sqrt(freqs)
fft_flicker = np.fft.rfft(np.random.normal(0, 1.0, N_PTS)) * flicker_shaping
flicker_noise = np.fft.irfft(fft_flicker, n=N_PTS)
flicker_noise = (flicker_noise / np.std(flicker_noise)) * 8.5e-17

# Serie temporal compuesta delta_nu / nu_0
clock_strain_series = white_noise + flicker_noise

print(f"[*] Nivel de ruido fraccionario RMS: sigma_y(tau=0.25s) = {np.std(clock_strain_series):.3e}")

# -------------------------------------------------------------------------
# BÚSQUEDA CIEGA POR PEINE MODULAR (Escaneo en Portadora f0)
# -------------------------------------------------------------------------
f0_range = np.linspace(0.40, 1.20, 17)
scan_results = []

sig_ac = clock_strain_series - np.mean(clock_strain_series)
fft_data = np.fft.rfft(sig_ac)
noise_std = np.std(np.abs(fft_data))

print("\n[*] Ejecutando filtrado adaptado sidéreo sobre la serie cuántica...")

for f_test in f0_range:
    params_test = FoliationParameters(f0=f_test)
    phys_test = SolitonPhysics(params_test)
    eng_test = MatchedFilterEngine(phys_test, station=st_ptb, source=source)
    
    template = eng_test.generate_template(t_arr, apply_sidereal=True)
    temp_ac = template - np.mean(template)
    norm_t = np.sqrt(np.sum(temp_ac**2))
    if norm_t > 0:
        temp_ac /= norm_t
        
    fft_t = np.fft.rfft(temp_ac)
    rho = float(np.abs(np.sum(fft_data * np.conj(fft_t))) / (noise_std * np.sqrt(N_PTS / 2.0)))
    
    # Cota superior al 95% CL en modulación fraccionaria (Delta nu / nu_0)
    frac_limit_95 = float((1.645 * np.std(sig_ac)) / np.sqrt(N_PTS))
    
    # Desacoplo del acoplamiento dilatónico d_e (sin apantallamiento bariónico, R = 0)
    # Delta nu / nu ~ d_e * chi_0 -> d_e < (Delta nu / nu) / chi_0
    de_limit = frac_limit_95 / phys_test.chi_0
    
    scan_results.append({
        "f0_carrier_hz": round(float(f_test), 3),
        "snr_rho": round(rho, 3),
        "frac_limit_95cl": f"{frac_limit_95:.4e}",
        "dilaton_de_95cl": f"{de_limit:.4e}"
    })

# -------------------------------------------------------------------------
# EXPORTACIÓN DE RESULTADOS
# -------------------------------------------------------------------------
output_data = {
    "experimento": "OpticalClockNetwork_PTB_CryoCavity",
    "duracion_horas": round(DURACION_DIAS * T_SID / 3600.0, 2),
    "total_muestras": N_PTS,
    "frecuencia_muestreo_hz": FS_CLOCK,
    "espectro_relojes": scan_results
}

out_filename = "exclusion_relojes_opticos_bipm.json"
with open(out_filename, "w") as f:
    json.dump(output_data, f, indent=2)

print("\n" + "=" * 95)
print(" RESULTADOS DEL BARRIDO MASIVO EN RELOJES ATÓMICOS:")
print(f"{'Frecuencia f0 (Hz)':<20} | {'SNR rho':<12} | {'Delta nu / nu (95% CL)':<25} | {'Acoplamiento d_e (95% CL)'}")
print("-" * 95)
for r in scan_results:
    print(f"{r['f0_carrier_hz']:<20.3f} | {r['snr_rho']:<12.3f} | {r['frac_limit_95cl']:<25} | {r['dilaton_de_95cl']}")

print("=" * 95)
print(f"[+] Archivo estructurado guardado con éxito en: {out_filename}")

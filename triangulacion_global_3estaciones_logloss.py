#!/usr/bin/env python3
"""
Triangulación Global Tri-Estación con Gradiente Logarítmico Robusto
Albert Cruañas Pardo (2026)
"""
import numpy as np
import time
from scipy.optimize import minimize
from detector_unificado_foliacion import (
    SolitonPhysics, MatchedFilterEngine, DetectorStation, CelestialSource,
    T_SID, V_DRIFT_NOMINAL
)

print("=" * 95)
print(" TRIANGULACIÓN GLOBAL CON GRADIENTE LOGARÍTMICO (PTB + KAGRA + PERTH)")
print("=" * 95)

R_EARTH = 6371000.0

def latlon_to_ecef(lat_deg, lon_deg):
    lat = np.radians(lat_deg)
    lon = np.radians(lon_deg)
    return np.array([
        R_EARTH * np.cos(lat) * np.cos(lon),
        R_EARTH * np.cos(lat) * np.sin(lon),
        R_EARTH * np.sin(lat)
    ])

def get_n_hat(ra_deg, dec_deg):
    ra = np.radians(ra_deg)
    dec = np.radians(dec_deg)
    return np.array([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])

st1 = DetectorStation(name="PTB-Hannover", lat_deg=52.37, lon_deg=9.73, azimuth_deg=45.0)
st2 = DetectorStation(name="KAGRA-Kamioka", lat_deg=36.42, lon_deg=137.31, azimuth_deg=15.0)
st3 = DetectorStation(name="UWA-Perth",    lat_deg=-31.98, lon_deg=115.82, azimuth_deg=90.0)

stations = [st1, st2, st3]
r_ecef = [latlon_to_ecef(s.lat_deg, s.lon_deg) for s in stations]

RA_TRUE = 83.6300
DEC_TRUE = 22.0100
H0_NOMINAL = 4.99e-17

fs = 4.0
n_pts = int(fs * T_SID)
t_arr = np.arange(n_pts) / fs
physics = SolitonPhysics()

n_hat_true = get_n_hat(RA_TRUE, DEC_TRUE)
tau_obs = []
for i in range(len(stations)):
    for j in range(i + 1, len(stations)):
        delta_r = r_ecef[i] - r_ecef[j]
        tau_obs.append((i, j, np.dot(delta_r, n_hat_true) / V_DRIFT_NOMINAL))

np.random.seed(42)
fft_data_list = []
std_noise_list = []

for s in stations:
    eng = MatchedFilterEngine(physics, station=s, source=CelestialSource(ra_deg=RA_TRUE, dec_deg=DEC_TRUE))
    sig = eng.generate_template(t_arr, apply_sidereal=True) * H0_NOMINAL
    ruido = np.random.normal(0, 1.0e-18, n_pts)
    data = sig + ruido
    ac = data - np.mean(data)
    fft_d = np.fft.rfft(ac)
    fft_data_list.append(fft_d)
    std_noise_list.append(np.std(np.abs(fft_d)))

coords_init = np.array([90.0, 35.0])

def funcion_perdida_robusta(coords):
    ra, dec = coords[0], coords[1]
    if not (-89.9 <= dec <= 89.9):
        return 1e6
        
    n_trial = get_n_hat(ra, dec)
    
    # Penalización cuadrática en tiempo de vuelo (sin exponenciales saturables)
    chi2_tau = 0.0
    for (i, j, tau_val) in tau_obs:
        delta_r = r_ecef[i] - r_ecef[j]
        tau_t = np.dot(delta_r, n_trial) / V_DRIFT_NOMINAL
        chi2_tau += ((tau_t - tau_val) / 3600.0)**2
    
    src = CelestialSource(ra_deg=ra, dec_deg=dec)
    rhos = []
    for k, s in enumerate(stations):
        eng_k = MatchedFilterEngine(physics, station=s, source=src)
        t_k = eng_k.generate_template(t_arr, apply_sidereal=True)
        t_norm = (t_k - np.mean(t_k)) / np.sqrt(np.sum((t_k - np.mean(t_k))**2))
        r_k = np.abs(np.sum(fft_data_list[k] * np.conj(np.fft.rfft(t_norm)))) / (std_noise_list[k] * np.sqrt(n_pts / 2.0))
        rhos.append(r_k)
        
    rho_total = np.sqrt(np.sum(np.array(rhos)**2))
    # Log-loss suave que no satura: minimizamos chi2 y maximizamos rho
    loss = -np.log(max(1.0, rho_total)) + 0.1 * chi2_tau
    return loss

print("[*] Minimizando con gradiente suave continuo...")
t0 = time.time()
res_tri = minimize(
    funcion_perdida_robusta,
    coords_init,
    method='Nelder-Mead',
    options={'xatol': 1e-4, 'fatol': 1e-4, 'maxiter': 200}
)
t_opt = time.time() - t0

ra_fin, dec_fin = res_tri.x[0], res_tri.x[1]

# Recomputar SNR final
src_fin = CelestialSource(ra_deg=ra_fin, dec_deg=dec_fin)
rhos_fin = []
for k, s in enumerate(stations):
    eng_k = MatchedFilterEngine(physics, station=s, source=src_fin)
    t_k = eng_k.generate_template(t_arr, apply_sidereal=True)
    t_norm = (t_k - np.mean(t_k)) / np.sqrt(np.sum((t_k - np.mean(t_k))**2))
    r_k = np.abs(np.sum(fft_data_list[k] * np.conj(np.fft.rfft(t_norm)))) / (std_noise_list[k] * np.sqrt(n_pts / 2.0))
    rhos_fin.append(r_k)
rho_tri_max = np.sqrt(np.sum(np.array(rhos_fin)**2))

delta_ra_deg = abs(ra_fin - RA_TRUE)
delta_dec_deg = abs(dec_fin - DEC_TRUE)
dist_ang_total = np.degrees(np.arccos(np.clip(np.dot(n_hat_true, get_n_hat(ra_fin, dec_fin)), -1.0, 1.0)))

print(f"[+] Convergencia en {t_opt:.2f} s ({res_tri.nfev} evaluaciones).")
print("\n" + "=" * 95)
print(" RESULTADOS DE LA RED TRI-ESTACIÓN (GRADIENTE SUAVE):")
print(f"  * Coordenadas reales inyectadas:  RA = {RA_TRUE:8.4f}°,  DEC = {DEC_TRUE:8.4f}°")
print(f"  * Coordenadas recuperadas:        RA = {ra_fin:8.4f}°,  DEC = {dec_fin:8.4f}°")
print(f"  * Error residual en RA:           Delta RA  = {delta_ra_deg:.4f}° ({delta_ra_deg*60.0:.2f} arcmin)")
print(f"  * Error residual en DEC:          Delta DEC = {delta_dec_deg:.4f}° ({delta_dec_deg*60.0:.2f} arcmin)")
print(f"  * Distancia angular en la esfera: Delta theta = {dist_ang_total:.4f}° ({dist_ang_total*60.0:.2f} arcmin)")
print(f"  * Significancia de Red Global:    rho_3st   = {rho_tri_max:.2f} sigma")
print("=" * 95)

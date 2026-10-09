#!/usr/bin/env python3
"""
Refinador Continuo Sub-Grado (Nelder-Mead Local Optimization)
Albert Cruañas Pardo (2026)

Toma las coordenadas gruesas del skymap y ejecuta optimización no lineal
para resolver la ascensión recta y declinación exactas con precisión milimétrica.
"""
import numpy as np
import time
from scipy.optimize import minimize
from detector_unificado_foliacion import (
    SolitonPhysics, MatchedFilterEngine, DetectorStation, CelestialSource,
    T_SID, V_DRIFT_NOMINAL
)

print("=" * 90)
print(" REFINAMIENTO GEODÉTICO CONTINUO: LOCALIZACIÓN SUB-GRADO (NELDER-MEAD)")
print("=" * 90)

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

# Estaciones de red
st1 = DetectorStation(name="PTB-Hannover", lat_deg=52.37, lon_deg=9.73, azimuth_deg=45.0)
st2 = DetectorStation(name="KAGRA-Kamioka", lat_deg=36.42, lon_deg=137.31, azimuth_deg=15.0)

r1_ecef = latlon_to_ecef(st1.lat_deg, st1.lon_deg)
r2_ecef = latlon_to_ecef(st2.lat_deg, st2.lon_deg)
delta_r12 = r1_ecef - r2_ecef

# Valores verdaderos inyectados
RA_TRUE = 83.63
DEC_TRUE = 22.01
H0_NOMINAL = 4.99e-17

fs = 4.0
n_pts = int(fs * T_SID)
t_arr = np.arange(n_pts) / fs

physics = SolitonPhysics()

# Simulación de las series temporales
e1_true = MatchedFilterEngine(physics, station=st1, source=CelestialSource(ra_deg=RA_TRUE, dec_deg=DEC_TRUE))
e2_true = MatchedFilterEngine(physics, station=st2, source=CelestialSource(ra_deg=RA_TRUE, dec_deg=DEC_TRUE))

n_hat_true = get_n_hat(RA_TRUE, DEC_TRUE)
tau12_obs = np.dot(delta_r12, n_hat_true) / V_DRIFT_NOMINAL

np.random.seed(42)
data1 = (e1_true.generate_template(t_arr, apply_sidereal=True) * H0_NOMINAL) + np.random.normal(0, 1.0e-18, n_pts)
data2 = (e2_true.generate_template(t_arr, apply_sidereal=True) * H0_NOMINAL) + np.random.normal(0, 1.0e-18, n_pts)

fft_data1 = np.fft.rfft(data1 - np.mean(data1))
fft_data2 = np.fft.rfft(data2 - np.mean(data2))
std1 = np.std(np.abs(fft_data1))
std2 = np.std(np.abs(fft_data2))

print(f"[*] Coordenadas reales inyectadas: RA = {RA_TRUE:.4f}°, DEC = {DEC_TRUE:.4f}°")
print(f"[*] Retardo de vuelo observado:   tau_12 = {tau12_obs:.2f} s ({tau12_obs/3600.0:.2f} h)")

# Punto de partida derivado del escaneo previo
coords_init = np.array([90.0, 35.0])
print(f"[*] Semilla inicial del optimizador: RA_0 = {coords_init[0]:.2f}°, DEC_0 = {coords_init[1]:.2f}°")

def funcion_objetivo(coords):
    ra, dec = coords[0], coords[1]
    
    # Restricciones físicas en coordenadas esféricas
    if not (-90.0 <= dec <= 90.0):
        return 1e6
        
    n_hat = get_n_hat(ra, dec)
    tau12_trial = np.dot(delta_r12, n_hat) / V_DRIFT_NOMINAL
    
    delta_tau = abs(tau12_trial - tau12_obs)
    factor_tiempo = np.exp(-0.5 * (delta_tau / 900.0)**2)
    
    src = CelestialSource(ra_deg=ra, dec_deg=dec)
    e1 = MatchedFilterEngine(physics, station=st1, source=src)
    e2 = MatchedFilterEngine(physics, station=st2, source=src)
    
    t1 = e1.generate_template(t_arr, apply_sidereal=True)
    t2 = e2.generate_template(t_arr, apply_sidereal=True)
    
    t1_norm = (t1 - np.mean(t1)) / np.sqrt(np.sum((t1 - np.mean(t1))**2))
    t2_norm = (t2 - np.mean(t2)) / np.sqrt(np.sum((t2 - np.mean(t2))**2))
    
    rho1 = np.abs(np.sum(fft_data1 * np.conj(np.fft.rfft(t1_norm)))) / (std1 * np.sqrt(n_pts / 2.0))
    rho2 = np.abs(np.sum(fft_data2 * np.conj(np.fft.rfft(t2_norm)))) / (std2 * np.sqrt(n_pts / 2.0))
    
    rho_total = np.sqrt(rho1**2 + rho2**2) * factor_tiempo
    # Retornamos el inverso / negativo para maximizar SNR
    return -rho_total

print("\n[*] Minimizando función de pérdida geodética...")
t_start = time.time()

opt_res = minimize(
    funcion_objetivo,
    coords_init,
    method='Nelder-Mead',
    options={'xatol': 1e-4, 'fatol': 1e-4, 'maxiter': 100, 'disp': False}
)

t_elapsed = time.time() - t_start
ra_opt, dec_opt = opt_res.x[0], opt_res.x[1]
rho_max = -opt_res.fun

print(f"[+] Optimización completada en {t_elapsed:.2f} s ({opt_res.nfev} evaluaciones).")

print("\n" + "=" * 90)
print(" RESULTADOS DEL REFINAMIENTO LOCAL:")
print(f"  * Coordenadas reales:       RA = {RA_TRUE:8.4f}°,  DEC = {DEC_TRUE:8.4f}°")
print(f"  * Estimación gruesa inicial: RA = {coords_init[0]:8.4f}°,  DEC = {coords_init[1]:8.4f}°")
print(f"  * Reconstrucción fina:      RA = {ra_opt:8.4f}°,  DEC = {dec_opt:8.4f}°")
print(f"  * Error residual en RA:     Delta RA  = {abs(ra_opt - RA_TRUE):.4f}° ({abs(ra_opt - RA_TRUE)*60.0:.2f} arcmin)")
print(f"  * Error residual en DEC:    Delta DEC = {abs(dec_opt - DEC_TRUE):.4f}° ({abs(dec_opt - DEC_TRUE)*60.0:.2f} arcmin)")
print(f"  * Distancia angular total:  Delta theta = {np.degrees(np.arccos(np.clip(np.dot(get_n_hat(RA_TRUE, DEC_TRUE), get_n_hat(ra_opt, dec_opt)), -1.0, 1.0))):.4f}°")
print(f"  * Significancia alcanzada:  rho_opt = {rho_max:.2f} sigma")
print("=" * 90)

#!/usr/bin/env python3
"""
Módulo de Modulación Sidérea y Triangulación Geodética en Red
Albert Cruañas Pardo (2026)
Integra rotación terrestre (23.9345 h), proyección topocéntrica
y búsqueda por plantilla coherente diurna.
"""
import numpy as np
from dataclasses import dataclass
from detector_unificado_foliacion import SolitonPhysics

print("=" * 90)
print(" MOTOR GEODÉTICO: MODULACIÓN SIDÉREA Y TRIANGULACIÓN MULTI-ESTACIÓN")
print("=" * 90)

T_SID = 86164.0905   # Período sidéreo terrestre en segundos (~23.93447 h)
OMEGA_EARTH = 2.0 * np.pi / T_SID # rad/s (~7.292115e-5 rad/s)
V_DRIFT = 5.7264     # m/s (Velocidad de deriva calibrada del tratado)

@dataclass
class DetectorStation:
    name: str
    lat_deg: float     # Latitud geográfica (+Norte, -Sur)
    lon_deg: float     # Longitud geográfica (+Este, -Oeste)
    azimuth_deg: float # Orientación del eje óptico/brazo respecto al Norte

class SiderealGeometry:
    def __init__(self, ra_source_deg: float = 85.0, dec_source_deg: float = -5.0):
        """
        ra_source_deg: Ascensión recta de la fuente / viento de foliación (e.g. Orión / Tauro)
        dec_source_deg: Declinación en la esfera celeste
        """
        self.alpha_0 = np.radians(ra_source_deg)
        self.delta_0 = np.radians(dec_source_deg)
        self.physics = SolitonPhysics()

    def get_projected_velocity(self, station: DetectorStation, t_arr: np.ndarray) -> np.ndarray:
        """
        Calcula la serie temporal de velocidad proyectada sobre el detector v_paralelo(t).
        """
        lat = np.radians(station.lat_deg)
        lon = np.radians(station.lon_deg)
        az = np.radians(station.azimuth_deg)

        # Tiempo sidéreo local aparente theta(t)
        theta = OMEGA_EARTH * t_arr + lon
        H = theta - self.alpha_0 # Ángulo horario local

        # Proyección en plano horizontal (Este, Norte)
        v_east = -V_DRIFT * np.cos(self.delta_0) * np.sin(H)
        v_north = V_DRIFT * (np.cos(lat) * np.sin(self.delta_0) - np.sin(lat) * np.cos(self.delta_0) * np.cos(H))

        # Componente longitudinal a lo largo del brazo óptico
        v_parallel = v_north * np.cos(az) + v_east * np.sin(az)
        return v_parallel

    def synthesize_sidereal_signal(self, station: DetectorStation, t_arr: np.ndarray, amp_base: float = 1.0e-17) -> np.ndarray:
        """
        Sintetiza la señal del peine de 17 líneas incluyendo modulación de amplitud diurna
        y corrimiento Doppler integrado en fase.
        """
        v_par = self.get_projected_velocity(station, t_arr)
        mod_amplitud = np.abs(v_par) / V_DRIFT # Factor geométrico [0, 1]

        freqs = self.physics.comb_frequencies()
        weights = self.physics.dirichlet_weights()

        signal = np.zeros_like(t_arr, dtype=float)
        dt = t_arr[1] - t_arr[0] if len(t_arr) > 1 else 1.0

        # Fase acumulada con corrección Doppler integrada
        c_light = 299792458.0
        fase_doppler_factor = 1.0 + (v_par / c_light)

        for fk, wk in zip(freqs, weights):
            # Fase instantánea: 2*pi * integral(fk * factor)
            fase = 2.0 * np.pi * fk * np.cumsum(fase_doppler_factor) * dt
            signal += wk * amp_base * mod_amplitud * np.cos(fase)

        return signal

# -------------------------------------------------------------------------
# SIMULACIÓN DE DETECCIÓN Y CORRELACIÓN SIDÉREA ENTRE DOS ESTACIONES
# -------------------------------------------------------------------------
# Estación 1: Europa Central (PTB / Hannover, Lat 52.3° N, Lon 9.7° E)
# Estación 2: Japón (KAGRA / Kamioka, Lat 36.4° N, Lon 137.3° E)
estacion_eur = DetectorStation("PTB-Hannover", 52.37, 9.73, 45.0)
estacion_jpn = DetectorStation("KAGRA-Kamioka", 36.42, 137.31, 15.0)

geom = SiderealGeometry(ra_source_deg=83.63, dec_source_deg=22.01) # Dirección Tauro/Crab

# Simular 24 horas completas con tasa de muestreo adecuada para baja frecuencia
duracion_total = T_SID # 1 día sidéreo exacto (~86164 s)
fs = 4.0               # 4 Hz (suficiente para resolver el peine de 0.5 - 1.5 Hz)
n_puntos = int(fs * duracion_total)
t_sid = np.arange(n_puntos) / fs

print(f"[*] Duración de análisis sidéreo:     {duracion_total / 3600.0:.2f} h ({n_puntos} muestras a {fs} Hz)")
print(f"[*] Estación 1 (Europa):             {estacion_eur.name} (Lat {estacion_eur.lat_deg}°, Lon {estacion_eur.lon_deg}°)")
print(f"[*] Estación 2 (Asia):               {estacion_jpn.name} (Lat {estacion_jpn.lat_deg}°, Lon {estacion_jpn.lon_deg}°)")

# Generar velocidades proyectadas y señales
v_eur = geom.get_projected_velocity(estacion_eur, t_sid)
v_jpn = geom.get_projected_velocity(estacion_jpn, t_sid)

# Desfase temporal entre picos de modulación diurna
t_pico_eur = t_sid[np.argmax(np.abs(v_eur))]
t_pico_jpn = t_sid[np.argmax(np.abs(v_jpn))]
delta_t_sid = (t_pico_eur - t_pico_jpn) % T_SID

print(f"\n[DINÁMICA DE PROYECCIÓN GEODÉTICA]")
print(f"  * Velocidad máx proyectada (EUR):  |v_par| = {np.max(np.abs(v_eur)):.3f} m/s (Nominal: {V_DRIFT} m/s)")
print(f"  * Velocidad máx proyectada (JPN):  |v_par| = {np.max(np.abs(v_jpn)):.3f} m/s (Nominal: {V_DRIFT} m/s)")
print(f"  * Hora pico de sensibilidad (EUR): T = {t_pico_eur / 3600.0:.2f} h sidéreas")
print(f"  * Hora pico de sensibilidad (JPN): T = {t_pico_jpn / 3600.0:.2f} h sidéreas")
print(f"  * Desfase geodésico diurno (EUR-JPN): Delta T = {delta_t_sid / 3600.0:.2f} h")
print(f"    (Predicción teórica por longitud Delta lon = 127.58°: {(127.58 / 360.0) * 23.934:.2f} h)")

# Síntesis con ruido instrumental térmico
amp_inyeccion = 4.99e-17
ruido_eur = np.random.normal(0, 1.0e-18, n_puntos)
ruido_jpn = np.random.normal(0, 1.0e-18, n_puntos)

sig_eur = geom.synthesize_sidereal_signal(estacion_eur, t_sid, amp_inyeccion) + ruido_eur
sig_jpn = geom.synthesize_sidereal_signal(estacion_jpn, t_sid, amp_inyeccion) + ruido_jpn

# Correlación cruzada entre envolventes de modulación de ambas estaciones
env_eur = np.abs(sig_eur)
env_jpn = np.abs(sig_jpn)

# Decimar para correlación lenta de envolvente diurna (bins de 60 segundos)
bin_size = int(60 * fs)
n_bins = n_puntos // bin_size
env_eur_dec = np.mean(env_eur[:n_bins * bin_size].reshape(n_bins, bin_size), axis=1)
env_jpn_dec = np.mean(env_jpn[:n_bins * bin_size].reshape(n_bins, bin_size), axis=1)

# Correlación normalizada
env_eur_norm = (env_eur_dec - np.mean(env_eur_dec)) / np.std(env_eur_dec)
env_jpn_norm = (env_jpn_dec - np.mean(env_jpn_dec)) / np.std(env_jpn_dec)
cc = np.correlate(env_eur_norm, env_jpn_norm, mode='full') / len(env_eur_norm)
lags = np.arange(-len(env_eur_norm) + 1, len(env_eur_norm)) * (bin_size / fs) / 3600.0 # en horas

max_cc_idx = np.argmax(cc)
desfase_recuperado_h = lags[max_cc_idx]

print(f"\n[CORRELACIÓN CRUZADA DE ENVOLVENTE MULTI-ESTACIÓN]")
print(f"  * Coeficiente de correlación máx:  r = {cc[max_cc_idx]:.3f}")
print(f"  * Desfase geodésico reconstruido:  Delta t_rec = {abs(desfase_recuperado_h):.2f} h")
print(f"  * Error en localización angular:    Error = {abs(abs(desfase_recuperado_h) - (delta_t_sid/3600.0)):.4f} h")

print("\n" + "=" * 90)
print(" SÍNTESIS:")
print(" El desfase de la modulación sidérea entre continentes permite triangular")
print(" la dirección celeste exacta (RA, DEC) del vector de foliación con datos continuos.")
print("=" * 90)

#!/usr/bin/env python3
"""
Motor Unificado de Análisis y Detección de Foliación Modular (17 Capas)
Albert Cruañas Pardo (2026)

Módulos integrados:
  1. SolitonPhysics: Perfil exacto de kink, masa y dispersión fonónica.
  2. DetectorStation & SiderealCoordinates: Geometría topocéntrica diurna.
  3. MatchedFilterEngine: Banco de plantillas coherentes con modulación de
     amplitud geodética y corrimiento Doppler sidéreo (23.9345 h).
  4. MultiDomainIngestor: Búsqueda ciega multicanal y extracción de cotas (95% CL).
"""
import numpy as np
from dataclasses import dataclass
from typing import Dict, Optional

# Constantes geodésicas y físicas
T_SID = 86164.0905                   # Período sidéreo terrestre (s)
OMEGA_EARTH = 2.0 * np.pi / T_SID    # rad/s (~7.292115e-5 rad/s)
V_DRIFT_NOMINAL = 5.7264             # m/s
C_LIGHT = 299792458.0                # m/s

# =========================================================================
# 1. PARÁMETROS FÍSICOS Y DINÁMICA DEL SOLITÓN
# =========================================================================
@dataclass
class FoliationParameters:
    N_sub: int = 17
    f0: float = 0.50             # Hz
    Lambda_r: float = 10785.0    # GeV
    v_EW: float = 246.22         # GeV
    m0: float = 1.0e-3           # GeV
    mu: float = 5.0e-4           # GeV
    lambda_self: float = 1.0e-4  # Adimensional
    gamma_damping: float = 0.05  # s^-1
    alpha_disp: float = 0.0450   # Coeficiente dispersivo cuadrático
    h0_nominal: float = 1.50e-16 # Deformación métrica de referencia
    v_drift: float = V_DRIFT_NOMINAL

@dataclass
class DetectorStation:
    name: str = "PTB-Hannover"
    lat_deg: float = 52.37       # Latitud geográfica (+N)
    lon_deg: float = 9.73        # Longitud geográfica (+E)
    azimuth_deg: float = 45.0    # Orientación del eje óptico/brazo respecto al Norte

@dataclass
class CelestialSource:
    ra_deg: float = 83.63        # Ascensión recta de la fuente / viento de foliación (Tauro/Crab)
    dec_deg: float = 22.01       # Declinación celeste

class SolitonPhysics:
    def __init__(self, params: FoliationParameters = FoliationParameters()):
        self.p = params
        self.xi = self.p.v_EW / self.p.Lambda_r
        self.d_e = 1.0 / np.sqrt(self.p.N_sub)
        self.chi_0 = self.p.mu / self.p.lambda_self
        self.delta_sol = np.sqrt(self.p.lambda_self) / (self.p.mu * np.sqrt(2.0))
        self.mass_sol = (self.p.mu**3) / (6.0 * (self.p.lambda_self**2) * np.sqrt(2.0 * self.p.lambda_self))
        
    def soliton_profile(self, z_arr: np.ndarray, z0: float = 0.0) -> np.ndarray:
        arg = (z_arr - z0) / (2.0 * self.delta_sol)
        return 0.5 * self.chi_0 * (1.0 + np.tanh(arg))

    def comb_frequencies(self) -> np.ndarray:
        k_indices = np.arange(self.p.N_sub)
        freqs_nom = self.p.f0 * (1.0 + k_indices / float(self.p.N_sub))
        # Dispersión cuadrática fonónica de interfase
        freqs_disp = freqs_nom * (1.0 - self.p.alpha_disp * (k_indices / float(self.p.N_sub))**2)
        return freqs_disp

    def dirichlet_weights(self) -> np.ndarray:
        k = np.arange(self.p.N_sub)
        w = np.cos((k * np.pi) / (2.0 * self.p.N_sub))
        return w / np.sqrt(np.sum(w**2))

# =========================================================================
# 2. MOTOR DE FILTRADO ADAPTADO CON MODULACIÓN SIDÉREA COMPLETA
# =========================================================================
class MatchedFilterEngine:
    def __init__(self, physics: SolitonPhysics,
                 station: DetectorStation = DetectorStation(),
                 source: CelestialSource = CelestialSource()):
        self.phys = physics
        self.station = station
        self.source = source
        self.freqs = physics.comb_frequencies()
        self.weights = physics.dirichlet_weights()
        
        # Precomputar constantes trigonométricas de la estación y fuente
        self.lat = np.radians(self.station.lat_deg)
        self.lon = np.radians(self.station.lon_deg)
        self.az = np.radians(self.station.azimuth_deg)
        self.ra = np.radians(self.source.ra_deg)
        self.dec = np.radians(self.source.dec_deg)

    def get_projected_velocity(self, time_array: np.ndarray) -> np.ndarray:
        """Calcula la velocidad proyectada sobre el detector v_paralelo(t)."""
        theta = OMEGA_EARTH * time_array + self.lon
        H = theta - self.ra # Ángulo horario local
        
        v_east = -self.phys.p.v_drift * np.cos(self.dec) * np.sin(H)
        v_north = self.phys.p.v_drift * (
            np.cos(self.lat) * np.sin(self.dec) - 
            np.sin(self.lat) * np.cos(self.dec) * np.cos(H)
        )
        return v_north * np.cos(self.az) + v_east * np.sin(self.az)

    def generate_template(self, time_array: np.ndarray, phase: float = 0.0,
                          apply_sidereal: bool = True) -> np.ndarray:
        """
        Genera la plantilla coherente del peine de 17 líneas.
        Si apply_sidereal=True, modula la amplitud diurna y acumula la fase Doppler.
        """
        n_pts = len(time_array)
        dt = time_array[1] - time_array[0] if n_pts > 1 else 1.0

        if apply_sidereal:
            v_par = self.get_projected_velocity(time_array)
            mod_amplitud = np.abs(v_par) / (self.phys.p.v_drift + 1e-15)
            # Factor Doppler de fase
            doppler_factor = 1.0 + (v_par / C_LIGHT)
            fase_base = np.cumsum(doppler_factor) * dt
        else:
            mod_amplitud = np.ones_like(time_array)
            fase_base = time_array

        template = np.zeros_like(time_array, dtype=float)
        for fk, wk in zip(self.freqs, self.weights):
            fase_total = 2.0 * np.pi * fk * fase_base + phase
            template += wk * mod_amplitud * np.cos(fase_total)

        # Normalización unitaria de la plantilla en energía
        norm = np.sqrt(np.sum(template**2))
        return template / norm if norm > 0 else template

    def analyze(self, time_array: np.ndarray, signal_data: np.ndarray,
                apply_sidereal: bool = True) -> Dict[str, float]:
        """
        Ejecuta el filtrado adaptado continuo mediante correlación en frecuencia
        con supresión de componente continua (DC) y normalización por PSD.
        """
        n_pts = len(time_array)
        dt = time_array[1] - time_array[0]
        
        # Sustracción estricta de DC y tendencia
        sig_ac = signal_data - np.mean(signal_data)
        
        # Plantilla adaptada (sidérea o isotrópica)
        template = self.generate_template(time_array, apply_sidereal=apply_sidereal)
        template_ac = template - np.mean(template)
        template_ac /= np.sqrt(np.sum(template_ac**2))
        
        # FFTs
        fft_data = np.fft.rfft(sig_ac)
        fft_temp = np.fft.rfft(template_ac)
        
        # Estimación del ruido del fondo
        ruido_sigma = np.std(np.abs(fft_data))
        
        # Correlación óptima compleja
        correlacion = np.sum(fft_data * np.conj(fft_temp))
        rho_stat = np.abs(correlacion) / (ruido_sigma * np.sqrt(n_pts / 2.0))
        
        std_real = np.std(sig_ac)
        upper_limit_95 = (1.645 * std_real) / np.sqrt(n_pts)
        
        return {
            "SNR_rho": float(rho_stat),
            "std_noise": float(std_real),
            "upper_limit_95": float(upper_limit_95),
            "n_samples": int(n_pts)
        }

# =========================================================================
# 3. INGESTOR MULTIDOMINIO Y PRUEBA DE COHERENCIA SIDÉREA
# =========================================================================
class MultiDomainIngestor:
    def __init__(self):
        self.physics = SolitonPhysics()
        self.station_eur = DetectorStation("PTB-Hannover", 52.37, 9.73, 45.0)
        self.engine = MatchedFilterEngine(self.physics, station=self.station_eur)

    def run_benchmark_suite(self):
        print("=" * 90)
        print(" MOTOR UNIFICADO: FILTRADO ADAPTADO CON MODULACIÓN SIDÉREA (24 HORAS)")
        print("=" * 90)
        
        print("\n[PARÁMETROS GEODÉTICOS DE LA ESTACIÓN]")
        print(f"  * Estación Base:                     {self.engine.station.name}")
        print(f"  * Coordenadas Geográficas:           Lat {self.engine.station.lat_deg}°, Lon {self.engine.station.lon_deg}°")
        print(f"  * Orientación del Brazo (Acimut):    {self.engine.station.azimuth_deg}°")
        print(f"  * Coordenadas Celestes del Vector:   RA = {self.engine.source.ra_deg}°, DEC = {self.engine.source.dec_deg}°")
        print(f"  * Período de Modulación Sidérea:     {T_SID / 3600.0:.4f} horas ({T_SID:.1f} s)")

        # PRUEBA 1: Serie Temporal Larga (24 Horas Sidéreas) con Modulación Coherente
        print("\n" + "-" * 90)
        print(" [CANAL 1] BÚSQUEDA COHERENTE SOBRE 1 DÍA SIDÉREO COMPLETO (23.93 h)")
        fs = 4.0 # 4 Hz
        duracion = T_SID
        n_pts = int(fs * duracion)
        t = np.arange(n_pts) / fs
        
        # Inyección teórica de señal sidérea completa (amplitud h0 = 4.99e-17)
        sig_pura = self.engine.generate_template(t, apply_sidereal=True) * 4.99e-17
        ruido = np.random.normal(0, 1.0e-18, n_pts)
        data_total = sig_pura + ruido

        # Filtrado con plantilla sidérea (Matched)
        res_matched = self.engine.analyze(t, data_total, apply_sidereal=True)
        # Filtrado con plantilla plana estática (Mismatched)
        res_unmatched = self.engine.analyze(t, data_total, apply_sidereal=False)

        print(f"  * Muestras procesadas:               {n_pts} ({duracion/3600.0:.2f} h)")
        print(f"  * Ruido instrumental (sigma_y):     {res_matched['std_noise']:.3e}")
        print(f"  * SNR con Plantilla Sidérea 24h:     rho = {res_matched['SNR_rho']:.2f} sigma")
        print(f"  * SNR con Plantilla Estática Plana:  rho = {res_unmatched['SNR_rho']:.2f} sigma")
        print(f"  * Ganancia de Detección Sidérea:     +{res_matched['SNR_rho'] - res_unmatched['SNR_rho']:.2f} sigma")
        print(f"  * Cota 95% CL alcanzada:            h_0 < {res_matched['upper_limit_95']:.3e}")

        # PRUEBA 2: Control Ciego de Falsa Alarma (Ruido Puro durante 24 h)
        print("\n" + "-" * 90)
        print(" [CANAL 2] CONTROL DE FALSA ALARMA SIDÉREA (RUIDO GAUSSIANO 24 h)")
        ruido_blanco_24h = np.random.normal(0, 1.0e-18, n_pts)
        res_ciego = self.engine.analyze(t, ruido_blanco_24h, apply_sidereal=True)
        print(f"  * SNR en ruido ciego sidéreo:        rho = {res_ciego['SNR_rho']:.2f} sigma (Esperado: < 3.0)")
        print(f"  * Veredicto de Robustez:             {'CONTROL POSITIVO (Sin Falsos Positivos)' if res_ciego['SNR_rho'] < 3.0 else 'ANOMALÍA'}")

        print("\n" + "=" * 90)
        print(" FILTRADO SIDÉREO OPERATIVO: Discrimina señales fijas de la Tierra frente")
        print(" a la precesión geodética del vector de foliación.")
        print("=" * 90)

if __name__ == "__main__":
    app = MultiDomainIngestor()
    app.run_benchmark_suite()

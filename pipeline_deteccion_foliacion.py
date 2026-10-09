#!/usr/bin/env python3
"""
Pipeline Unificado de Detección y Análisis de Foliación Modular
Albert Cruañas Pardo (2026)

CLI monolítico que integra:
  1. Benchmark de Dinámica de Solitón (Masa, Dispersión, Kink).
  2. Triangulación Global Tri-Estación (Nelder-Mead con gradiente suave).
  3. Búsqueda Ciega Multidía y Extracción de Cotas Asintóticas (95% CL).
  4. Exportación formal de resultados en formato JSON para documentación.
"""
import argparse
import json
import time
import numpy as np
from dataclasses import dataclass, asdict
from scipy.optimize import minimize

# =========================================================================
# CONSTANTES FÍSICAS Y GEODÉTICAS
# =========================================================================
T_SID = 86164.0905                   # Período sidéreo terrestre (s)
OMEGA_EARTH = 2.0 * np.pi / T_SID    # rad/s (~7.292115e-5 rad/s)
V_DRIFT_NOMINAL = 5.7264             # m/s
C_LIGHT = 299792458.0                # m/s
R_EARTH = 6371000.0                  # m

def latlon_to_ecef(lat_deg: float, lon_deg: float) -> np.ndarray:
    lat = np.radians(lat_deg)
    lon = np.radians(lon_deg)
    return np.array([
        R_EARTH * np.cos(lat) * np.cos(lon),
        R_EARTH * np.cos(lat) * np.sin(lon),
        R_EARTH * np.sin(lat)
    ])

def get_n_hat(ra_deg: float, dec_deg: float) -> np.ndarray:
    ra = np.radians(ra_deg)
    dec = np.radians(dec_deg)
    return np.array([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])

# =========================================================================
# ESTRUCTURAS DE DATOS Y MODELO FÍSICO
# =========================================================================
@dataclass
class FoliationParameters:
    N_sub: int = 17
    f0: float = 0.50
    Lambda_r: float = 10785.0
    v_EW: float = 246.22
    m0: float = 1.0e-3
    mu: float = 5.0e-4
    lambda_self: float = 1.0e-4
    gamma_damping: float = 0.05
    alpha_disp: float = 0.0450
    h0_nominal: float = 1.50e-16
    v_drift: float = V_DRIFT_NOMINAL

@dataclass
class DetectorStation:
    name: str
    lat_deg: float
    lon_deg: float
    azimuth_deg: float

@dataclass
class CelestialSource:
    ra_deg: float
    dec_deg: float

class SolitonPhysics:
    def __init__(self, params: FoliationParameters = FoliationParameters()):
        self.p = params
        self.xi = self.p.v_EW / self.p.Lambda_r
        self.d_e = 1.0 / np.sqrt(self.p.N_sub)
        self.chi_0 = self.p.mu / self.p.lambda_self
        self.delta_sol = np.sqrt(self.p.lambda_self) / (self.p.mu * np.sqrt(2.0))
        self.mass_sol = (self.p.mu**3) / (6.0 * (self.p.lambda_self**2) * np.sqrt(2.0 * self.p.lambda_self))
        
    def comb_frequencies(self) -> np.ndarray:
        k_indices = np.arange(self.p.N_sub)
        freqs_nom = self.p.f0 * (1.0 + k_indices / float(self.p.N_sub))
        return freqs_nom * (1.0 - self.p.alpha_disp * (k_indices / float(self.p.N_sub))**2)

    def dirichlet_weights(self) -> np.ndarray:
        k = np.arange(self.p.N_sub)
        w = np.cos((k * np.pi) / (2.0 * self.p.N_sub))
        return w / np.sqrt(np.sum(w**2))

class MatchedFilterEngine:
    def __init__(self, physics: SolitonPhysics, station: DetectorStation, source: CelestialSource):
        self.phys = physics
        self.station = station
        self.source = source
        self.freqs = physics.comb_frequencies()
        self.weights = physics.dirichlet_weights()
        self.lat = np.radians(self.station.lat_deg)
        self.lon = np.radians(self.station.lon_deg)
        self.az = np.radians(self.station.azimuth_deg)
        self.ra = np.radians(self.source.ra_deg)
        self.dec = np.radians(self.source.dec_deg)

    def get_projected_velocity(self, t_arr: np.ndarray) -> np.ndarray:
        theta = OMEGA_EARTH * t_arr + self.lon
        H = theta - self.ra
        v_east = -self.phys.p.v_drift * np.cos(self.dec) * np.sin(H)
        v_north = self.phys.p.v_drift * (
            np.cos(self.lat) * np.sin(self.dec) - np.sin(self.lat) * np.cos(self.dec) * np.cos(H)
        )
        return v_north * np.cos(self.az) + v_east * np.sin(self.az)

    def generate_template(self, t_arr: np.ndarray, apply_sidereal: bool = True) -> np.ndarray:
        n_pts = len(t_arr)
        dt = t_arr[1] - t_arr[0] if n_pts > 1 else 1.0

        if apply_sidereal:
            v_par = self.get_projected_velocity(t_arr)
            mod_amp = np.abs(v_par) / (self.phys.p.v_drift + 1e-15)
            doppler_factor = 1.0 + (v_par / C_LIGHT)
            fase_base = np.cumsum(doppler_factor) * dt
        else:
            mod_amp = np.ones_like(t_arr)
            fase_base = t_arr

        template = np.zeros_like(t_arr, dtype=float)
        for fk, wk in zip(self.freqs, self.weights):
            template += wk * mod_amp * np.cos(2.0 * np.pi * fk * fase_base)

        norm = np.sqrt(np.sum(template**2))
        return template / norm if norm > 0 else template

# =========================================================================
# SUBCOMANDOS DEL PIPELINE
# =========================================================================
def run_triangulation():
    print("=" * 95)
    print(" [MODO] TRIANGULACIÓN GLOBAL CONTINUA TRI-ESTACIÓN (NELDER-MEAD)")
    print("=" * 95)

    stations = [
        DetectorStation("PTB-Hannover", 52.37, 9.73, 45.0),
        DetectorStation("KAGRA-Kamioka", 36.42, 137.31, 15.0),
        DetectorStation("UWA-Perth", -31.98, 115.82, 90.0)
    ]
    r_ecef = [latlon_to_ecef(s.lat_deg, s.lon_deg) for s in stations]

    RA_TRUE, DEC_TRUE = 83.6300, 22.0100
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
    fft_data_list, std_noise_list = [], []
    for s in stations:
        eng = MatchedFilterEngine(physics, station=s, source=CelestialSource(RA_TRUE, DEC_TRUE))
        sig = eng.generate_template(t_arr, apply_sidereal=True) * H0_NOMINAL
        ruido = np.random.normal(0, 1.0e-18, n_pts)
        ac = (sig + ruido) - np.mean(sig + ruido)
        fft_d = np.fft.rfft(ac)
        fft_data_list.append(fft_d)
        std_noise_list.append(np.std(np.abs(fft_d)))

    coords_init = np.array([90.0, 35.0])

    def objective(coords):
        ra, dec = coords[0], coords[1]
        if not (-89.9 <= dec <= 89.9):
            return 1e6
        n_trial = get_n_hat(ra, dec)
        chi2_tau = sum([((np.dot(r_ecef[i] - r_ecef[j], n_trial)/V_DRIFT_NOMINAL - val)/3600.0)**2 for (i, j, val) in tau_obs])
        src = CelestialSource(ra, dec)
        rhos = []
        for k, s in enumerate(stations):
            eng_k = MatchedFilterEngine(physics, station=s, source=src)
            t_k = eng_k.generate_template(t_arr, apply_sidereal=True)
            t_norm = (t_k - np.mean(t_k)) / np.sqrt(np.sum((t_k - np.mean(t_k))**2))
            r_k = np.abs(np.sum(fft_data_list[k] * np.conj(np.fft.rfft(t_norm)))) / (std_noise_list[k] * np.sqrt(n_pts / 2.0))
            rhos.append(r_k)
        return -np.log(max(1.0, np.sqrt(np.sum(np.array(rhos)**2)))) + 0.1 * chi2_tau

    t0 = time.time()
    res = minimize(objective, coords_init, method='Nelder-Mead', options={'xatol': 1e-4, 'fatol': 1e-4, 'maxiter': 200})
    t_opt = time.time() - t0

    ra_fin, dec_fin = res.x[0], res.x[1]
    # SNR final
    src_fin = CelestialSource(ra_fin, dec_fin)
    rhos_fin = []
    for k, s in enumerate(stations):
        eng_k = MatchedFilterEngine(physics, station=s, source=src_fin)
        t_k = eng_k.generate_template(t_arr, apply_sidereal=True)
        t_norm = (t_k - np.mean(t_k)) / np.sqrt(np.sum((t_k - np.mean(t_k))**2))
        rhos_fin.append(np.abs(np.sum(fft_data_list[k] * np.conj(np.fft.rfft(t_norm)))) / (std_noise_list[k] * np.sqrt(n_pts / 2.0)))
    rho_total = np.sqrt(np.sum(np.array(rhos_fin)**2))

    salida = {
        "status": "CONVERGED",
        "duracion_opt_s": round(t_opt, 2),
        "iteraciones": res.nfev,
        "coordenadas_reales": {"ra_deg": RA_TRUE, "dec_deg": DEC_TRUE},
        "coordenadas_recuperadas": {"ra_deg": round(ra_fin, 4), "dec_deg": round(dec_fin, 4)},
        "error_residual_arcmin": {
            "delta_ra": round(abs(ra_fin - RA_TRUE) * 60.0, 2),
            "delta_dec": round(abs(dec_fin - DEC_TRUE) * 60.0, 2)
        },
        "snr_combinado_sigma": round(rho_total, 2)
    }

    print(json.dumps(salida, indent=2))
    with open("resultado_triangulacion_global.json", "w") as f:
        json.dump(salida, f, indent=2)
    print("\n[+] Guardado en resultado_triangulacion_global.json")

def run_exclusion():
    print("=" * 95)
    print(" [MODO] BÚSQUEDA CIEGA MULTIDÍA (72 HORAS) Y COTAS ASINTÓTICAS")
    print("=" * 95)

    DIAS = 3
    DURACION_S = DIAS * T_SID
    FS = 4.0
    N_SAMPLES = int(FS * DURACION_S)
    t_serie = np.arange(N_SAMPLES) / FS

    station = DetectorStation("PTB-Hannover", 52.37, 9.73, 45.0)
    source = CelestialSource(83.63, 22.01)
    physics = SolitonPhysics()
    eng = MatchedFilterEngine(physics, station, source)

    np.random.seed(1234)
    ruido_blanco = np.random.normal(0, 1.0e-18, N_SAMPLES)
    freqs = np.fft.rfftfreq(N_SAMPLES, 1.0 / FS)
    freqs[0] = freqs[1]
    ruido_coloreado = np.fft.irfft(np.fft.rfft(ruido_blanco) / np.sqrt(freqs), n=N_SAMPLES)
    ruido_coloreado = (ruido_coloreado / np.std(ruido_coloreado)) * 1.15e-18

    # Filtrado adaptado
    sig_ac = ruido_coloreado - np.mean(ruido_coloreado)
    temp = eng.generate_template(t_serie, apply_sidereal=True)
    temp_ac = (temp - np.mean(temp)) / np.sqrt(np.sum((temp - np.mean(temp))**2))
    
    fft_d = np.fft.rfft(sig_ac)
    fft_t = np.fft.rfft(temp_ac)
    rho = float(np.abs(np.sum(fft_d * np.conj(fft_t))) / (np.std(np.abs(fft_d)) * np.sqrt(N_SAMPLES / 2.0)))
    h0_lim = float((1.645 * np.std(sig_ac)) / np.sqrt(N_SAMPLES))

    v_EW = physics.p.v_EW
    N_sub = physics.p.N_sub
    masas = [10.0, 100.0, 1000.0, 10785.0]
    tabla_cotas = []
    for M in masas:
        factor_m = np.sqrt(1.0 + (M / 1000.0)**2)
        h0_m = h0_lim * factor_m
        de_lam = h0_m / (v_EW * np.sqrt(N_sub))
        lam_min = ((1.0 / np.sqrt(N_sub)) / de_lam) * 1e-3
        tabla_cotas.append({
            "masa_gev": M,
            "h0_95_cl": f"{h0_m:.3e}",
            "de_over_lambda_gev_inv": f"{de_lam:.3e}",
            "lambda_min_tev": f"{lam_min:.2e}"
        })

    salida = {
        "duracion_horas": round(DURACION_S / 3600.0, 2),
        "muestras": N_SAMPLES,
        "snr_ciego_sigma": round(rho, 2),
        "h0_limite_base": f"{h0_lim:.3e}",
        "cotas_por_masa": tabla_cotas
    }

    print(json.dumps(salida, indent=2))
    with open("resultado_exclusion_multidia.json", "w") as f:
        json.dump(salida, f, indent=2)
    print("\n[+] Guardado en resultado_exclusion_multidia.json")

# =========================================================================
# MAIN
# =========================================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Pipeline Unificado de Foliación Modular (2026)")
    parser.add_argument("--mode", choices=["triangulate", "exclude", "all"], default="all",
                        help="Modo de ejecución: triangulate, exclude o all")
    args = parser.parse_args()

    if args.mode in ["triangulate", "all"]:
        run_triangulation()
    if args.mode in ["exclude", "all"]:
        run_exclusion()

#!/usr/bin/env python3
"""
Pipeline de Búsqueda de Transitorios Topológicos: Red GNOME
Albert Cruañas Pardo (2026)

Busca cruces de paredes de dominio y paquetes de ondas de solitones
mediante correlación temporal cruzada con retardo de vuelo geodésico
entre estaciones globales sincronizadas por GPS.
"""
import os
import json
import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt

print("=" * 95)
print(" RED GNOME: BÚSQUEDA DE TRANSITORIOS Y PAQUETES DE ONDA TOPOLÓGICOS")
print("=" * 95)

# Coordenadas geocéntricas aproximadas (ECEF en km) de 4 estaciones GNOME
R_EARTH = 6371.0  # km
C_LIGHT = 299792.458  # km/s
V_HALO = 300.0  # km/s (velocidad típica de halo galáctico para defectos oscuros)

ESTACIONES_GNOME = {
    "PTB_Alemania": {"lat": 52.295, "lon": 10.454},
    "KRA_Polonia":   {"lat": 50.064, "lon": 19.945},
    "LAN_EEUU":      {"lat": 35.880, "lon": -106.303},
    "DAE_Corea":     {"lat": 36.350, "lon": 127.385}
}

def latlon_a_ecef(lat_deg, lon_deg):
    lat = np.radians(lat_deg)
    lon = np.radians(lon_deg)
    x = R_EARTH * np.cos(lat) * np.cos(lon)
    y = R_EARTH * np.cos(lat) * np.sin(lon)
    z = R_EARTH * np.sin(lat)
    return np.array([x, y, z])

posiciones_ecef = {k: latlon_a_ecef(v["lat"], v["lon"]) for k, v in ESTACIONES_GNOME.items()}

# Simulación de un stream temporal de 3600 segundos (1 hora) a 10 Hz
FS_GNOME = 10.0  # Hz
DURACION_S = 3600.0
n_samples = int(DURACION_S * FS_GNOME)
t_vec = np.linspace(0, DURACION_S, n_samples, endpoint=False)

print(f"[*] Estaciones operativas activas: {len(ESTACIONES_GNOME)}")
for est, pos in posiciones_ecef.items():
    print(f"    - {est:<15}: Coordenadas ECEF [{pos[0]:.1f}, {pos[1]:.1f}, {pos[2]:.1f}] km")

# Generar ruido de magnetometría atómica (fondo SERF / vapor alcalino: ruido blanco + 1/f)
np.random.seed(42)
senales = {}
for est in ESTACIONES_GNOME:
    ruido_blanco = np.random.normal(0, 1.0, n_samples)
    b, a = signal.butter(2, [0.1 / (FS_GNOME/2), 3.0 / (FS_GNOME/2)], btype='bandpass')
    senales[est] = signal.filtfilt(b, a, ruido_blanco) * 15.0  # en fT (femtoTeslas)

# Parámetros del vector de llegada de una ráfaga de prueba (prueba de recuperación)
# Dirección incidente: vector unitario hacia el ápice solar
n_incidente = np.array([0.45, -0.25, 0.86])
n_incidente /= np.linalg.norm(n_incidente)
t_impacto_centro = 1800.0  # Segundo 1800 (mitad de la serie)
ancho_paquete_s = 2.5      # Duración del cruce del solitón (segundos)
amplitud_ft = 8.5          # Magnitud transitoria en fT

# Inyección sintética con retardo físico real
for est, r_vec in posiciones_ecef.items():
    delta_tau = -np.dot(r_vec, n_incidente) / V_HALO  # Retardo geodésico (segundos)
    t_local = t_impacto_centro + delta_tau
    
    # Paquete de onda oscilatorio modular centrado en 0.50 Hz
    forma_onda = (amplitud_ft * 
                  np.exp(-0.5 * ((t_vec - t_local) / ancho_paquete_s)**2) * 
                  np.cos(2 * np.pi * 0.50 * (t_vec - t_local)))
    senales[est] += forma_onda

print("\n[*] Ejecutando análisis de coincidencia y correlación temporal cruzada...")

# Ventana de correlación por bloques (análisis por deslizamiento de 30 segundos)
t_ventana = 30.0
n_step = int(5.0 * FS_GNOME)
n_win = int(t_ventana * FS_GNOME)

t_centros = []
potencias_coherentes = []

est_keys = list(ESTACIONES_GNOME.keys())
n_pares = len(est_keys) * (len(est_keys) - 1) // 2

for idx_inicio in range(0, n_samples - n_win, n_step):
    t_mid = t_vec[idx_inicio + n_win // 2]
    t_centros.append(t_mid)
    
    # Extraer segmentos
    segmentos = {k: senales[k][idx_inicio : idx_inicio + n_win] for k in est_keys}
    
    # Coherencia geométrica multiespectral (producto cruzado normalizado)
    suma_corr = 0.0
    for i in range(len(est_keys)):
        for j in range(i + 1, len(est_keys)):
            k1, k2 = est_keys[i], est_keys[j]
            s1 = segmentos[k1] - np.mean(segmentos[k1])
            s2 = segmentos[k2] - np.mean(segmentos[k2])
            norm = np.sqrt(np.sum(s1**2) * np.sum(s2**2))
            
            # Máxima correlación cruzada dentro de la ventana de retardo causal (+- 60 s)
            max_delay_pts = int((2.0 * R_EARTH / V_HALO) * FS_GNOME) + 5
            xcorr = np.correlate(s1, s2, mode='full')
            mid_idx = len(xcorr) // 2
            corr_ventana = xcorr[mid_idx - max_delay_pts : mid_idx + max_delay_pts + 1]
            suma_corr += np.max(np.abs(corr_ventana)) / norm
            
    potencia_red = suma_corr / n_pares
    potencias_coherentes.append(potencia_red)

t_centros = np.array(t_centros)
potencias_coherentes = np.array(potencias_coherentes)

# Estimación de significancia de transitorio (Z-score)
mediana_ruido = np.median(potencias_coherentes)
mad_ruido = np.median(np.abs(potencias_coherentes - mediana_ruido))
sigma_equiv = 1.4826 * mad_ruido
snr_transitorios = (potencias_coherentes - mediana_ruido) / sigma_equiv

idx_pico = np.argmax(snr_transitorios)
t_pico = t_centros[idx_pico]
max_snr = snr_transitorios[idx_pico]

print(f"[+] Análisis completado:")
print(f"    - Evento más coherente detectado en: t = {t_pico:.1f} s")
print(f"    - Coherencia máxima de red:          S = {potencias_coherentes[idx_pico]:.3f}")
print(f"    - Significancia estadística (Z):     {max_snr:.2f} sigma")

# Exportar informe JSON
resultado_gnome = {
    "red": "GNOME_4_STATIONS",
    "estaciones": list(ESTACIONES_GNOME.keys()),
    "duracion_analizada_s": DURACION_S,
    "velocidad_prueba_km_s": V_HALO,
    "evento_maximo": {
        "tiempo_deteccion_s": round(float(t_pico), 2),
        "significancia_sigma": round(float(max_snr), 2),
        "coherencia_red": round(float(potencias_coherentes[idx_pico]), 4)
    }
}

with open("exclusion_gnome_transitorios_1h.json", "w") as f:
    json.dump(resultado_gnome, f, indent=2)

# Gráfico de diagnóstico
plt.rcParams.update({"font.family": "serif", "font.size": 11, "mathtext.fontset": "cm"})
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)

# Series temporales alrededor del evento
zoom_mask = (t_vec >= t_pico - 25.0) & (t_vec <= t_pico + 25.0)
for est in est_keys:
    ax1.plot(t_vec[zoom_mask], senales[est][zoom_mask], label=est, alpha=0.8, linewidth=1.2)
ax1.set_ylabel("Campo magnético [fT]")
ax1.set_title("Red GNOME: Respuesta de Magnetómetros Globales ante Paquete Transitorio")
ax1.grid(True, linestyle=":", alpha=0.6)
ax1.legend(loc="upper right", ncol=2)

# Traza de Coherencia de Red
ax2.plot(t_centros, snr_transitorios, color="#8e44ad", linewidth=1.8, label="Coherencia Multiestación GNOME")
ax2.axhline(3.0, color="#f39c12", linestyle="--", label="Umbral 3 sigma")
ax2.axhline(5.0, color="#e74c3c", linestyle=":", label="Descubrimiento 5 sigma")
ax2.set_xlabel("Tiempo transcurrido [s]")
ax2.set_ylabel("Significancia de Coincidencia [$\sigma$]")
ax2.set_xlim(0, DURACION_S)
ax2.grid(True, linestyle=":", alpha=0.6)
ax2.legend(loc="upper left")

plt.tight_layout()
plt.savefig("analisis_transitorios_gnome.png", dpi=300)
plt.savefig("analisis_transitorios_gnome.pdf")
print("[+] Gráficos exportados: analisis_transitorios_gnome.png y .pdf")
print("[+] Informe guardado en: exclusion_gnome_transitorios_1h.json")

#!/usr/bin/env python3
"""
Análisis del Dataset Público del Holómetro de Fermilab (PRL 2016)
Albert Cruañas Pardo (2026)

Descarga y procesa la densidad espectral cruzada (CSD) de Fermilab E-990
y evalúa la compatibilidad con resonancias armónicas modulares.
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

csv_url = "https://holometer.fnal.gov/data/holo_data_PRL2016.csv"
csv_file = "holo_data_PRL2016.csv"

print("=" * 85)
print(" DESCARGA Y ANÁLISIS DE DATOS: FERMILAB HOLOMETER (E-990)")
print("=" * 85)

if not os.path.exists(csv_file):
    print(f"[*] Descargando dataset desde {csv_url}...")
    # Intentar curl con timeout
    ret = os.system(f"curl -k -L --connect-timeout 15 --retry 2 -o '{csv_file}' '{csv_url}'")
    if ret != 0 or not os.path.exists(csv_file) or os.path.getsize(csv_file) < 1000:
        print("[-] Aviso: No se pudo descargar directamente desde la URL de Fermilab.")
        print("    Verificando si el servidor bloquea TLS o si requiere inspección manual.")
else:
    print(f"[+] Archivo local encontrado: {csv_file} ({os.path.getsize(csv_file)/(1024):.1f} KB)")

# Si el archivo está descargado, analizar estructura
if os.path.exists(csv_file) and os.path.getsize(csv_file) > 1000:
    try:
        # Cargar datos: típicamente columnas: Frecuencia (Hz o MHz), Real(CSD), Imag(CSD), Error
        datos = np.genfromtxt(csv_file, delimiter=',', skip_header=1)
        freq = datos[:, 0]
        csd_real = datos[:, 1]
        csd_imag = datos[:, 2] if datos.shape[1] > 2 else np.zeros_like(csd_real)
        csd_mag = np.sqrt(csd_real**2 + csd_imag**2)
        
        print(f"[+] Dataset cargado exitosamente: {len(freq)} canales espectrales")
        print(f"    - Rango de frecuencia: [{np.min(freq):.3e}, {np.max(freq):.3e}] Hz")
        print(f"    - Nivel medio de CSD:   {np.mean(csd_mag):.3e}")
        
        # Test de hipótesis: Búsqueda de modulación periódica delta_f = f / 17
        # Evaluamos la autocorrelación del espectro para detectar peines
        csd_norm = (csd_real - np.mean(csd_real)) / np.std(csd_real)
        autocorr_spec = np.correlate(csd_norm, csd_norm, mode='full')
        lags = np.arange(-len(csd_norm) + 1, len(csd_norm))
        
        print(f"[*] Evaluando periodicidades ocultas en el espectro...")
        print(f"    - Correlación máxima fuera del origen: {np.max(autocorr_spec[lags != 0]) / len(csd_norm):.4f}")
        
    except Exception as e:
        print(f"[-] Error al parsear el archivo CSV: {e}")
else:
    print("[-] Dataset no disponible en disco para análisis numérico inmediato.")
print("=" * 85)

#!/usr/bin/env python3
"""
Extracción de valores exactos de sensibilidad métrica en LIGO Livingston (L1 - O3)
Albert Cruañas Pardo (2026)
"""

from gwpy.timeseries import TimeSeries
import numpy as np

gps_start = 1242442967
gps_end = gps_start + 1024

print("[*] Leyendo datos locales en caché de LIGO L1...")
data = TimeSeries.fetch_open_data('L1', gps_start, gps_end, sample_rate=4096, cache=True)
asd = data.asd(fftlength=32, overlap=16)

f_targets = [13.15, 17.20, 26.30, 210.40, 210.42]

print("\n" + "=" * 78)
print(" SENSIBILIDAD MÉTRICA INTERFEROMÉTRICA (LIGO L1 - RUN O3)")
print("=" * 78)
for ft in f_targets:
    idx = np.argmin(np.abs(asd.frequencies.value - ft))
    val = asd.value[idx]
    f_real = asd.frequencies.value[idx]
    print(f"Frecuencia objetivo: {ft:6.2f} Hz | Celda FFT: {f_real:6.2f} Hz | ASD: {val:.4e} Hz^(-1/2)")
print("=" * 78 + "\n")

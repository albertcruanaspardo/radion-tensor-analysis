#!/usr/bin/env python3
"""
Análisis Espectral de Deformación Métrica con Datos Abiertos de LIGO (O3)
Albert Cruañas Pardo (2026)
"""

import matplotlib.pyplot as plt
from gwpy.timeseries import TimeSeries

# Intervalo GPS representativo del run O3 (1024 segundos de observación continua)
gps_start = 1242442967
gps_end = gps_start + 1024

print("=" * 80)
print(" DESCARGA Y ANÁLISIS DE DATOS INTERFEROMÉTRICOS (GWOSC / LIGO O3)")
print("=" * 80)
print(f"Detector         : LIGO Livingston (L1)")
print(f"Ventana GPS      : {gps_start} a {gps_end} (1024 segundos)")
print("=" * 80)

print("[*] Conectando con los servidores de GWOSC y descargando canal strain h(t)...")
data = TimeSeries.fetch_open_data('L1', gps_start, gps_end, sample_rate=4096, cache=True)
print("[+] Datos descargados con éxito.")

# Cálculo de la densidad espectral de amplitud (ASD)
print("[*] Calculando Densidad Espectral de Amplitud (ASD)...")
asd = data.asd(fftlength=32, overlap=16)

# Marcadores de frecuencia a contrastar
f_radion = 13.15
f_sonda  = 210.42
f_arm16  = 16 * 13.15  # 210.40 Hz

# Graficar
plt.figure(figsize=(11, 6))
plt.plot(asd.frequencies.value, asd.value, color='black', lw=0.7, label='LIGO L1 (Run O3)')

plt.axvline(f_radion, color='green', linestyle='--', lw=1.2, label=f'Bomba radión ({f_radion:.2f} Hz)')
plt.axvline(f_arm16, color='red', linestyle=':', lw=1.2, label=f'Armónico 16 ({f_arm16:.2f} Hz)')
plt.axvline(f_sonda, color='cyan', linestyle='-.', lw=1.2, label=f'Sonda portadora ({f_sonda:.2f} Hz)')

plt.xlim(10, 250)
plt.yscale('log')
plt.xlabel('Frecuencia [Hz]', fontsize=10)
plt.ylabel(r'Sensibilidad de Deformación $[1/\sqrt{\text{Hz}}]$', fontsize=10)
plt.title('Espectro de Ruido y Deformación Métrica: LIGO Livingston (O3)', fontweight='bold')
plt.grid(True, which='both', linestyle=':', alpha=0.6)
plt.legend(loc='upper right')

out_png = "ligo_espectro_test.png"
plt.tight_layout()
plt.savefig(out_png, dpi=300)
print(f"[+] Gráfica generada y guardada en: {out_png}\n")

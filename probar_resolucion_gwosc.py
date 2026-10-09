#!/usr/bin/env python3
"""
Test de Resolución Espectral y Enventanado Anti-Fuga para GWOSC
Albert Cruañas Pardo (2026)
"""
import numpy as np
from scipy.signal import windows
from detector_unificado_foliacion import SolitonPhysics, MatchedFilterEngine

print("=" * 85)
print(" DIAGNÓSTICO DE RESOLUCIÓN ESPECTRAL Y CONTROL DE FUGA (LEAKAGE)")
print("=" * 85)

physics = SolitonPhysics()
engine = MatchedFilterEngine(physics)

# Comparativa entre dos duraciones: 32 s (baja resolución) vs 256 s (resolución adecuada)
fs = 16.0
delta_f_peine = 0.50 / 17.0 # 0.0294 Hz

for T in [32.0, 256.0]:
    N = int(fs * T)
    t = np.arange(N) / fs
    df = 1.0 / T
    
    # Ruido blanco puro pasado por la misma banda
    ruido = np.random.normal(0, 1e-21, N)
    
    # Aplicar ventana Tukey (tapering del 10% en los extremos para eliminar discontinuidad)
    win = windows.tukey(N, alpha=0.1)
    ruido_enventanado = ruido * win
    
    res = engine.analyze(t, ruido_enventanado)
    
    print(f"\n[Duración T = {T:.0f} s | Muestras = {N}]")
    print(f"  * Resolución espectral (Delta f):    {df:.5f} Hz")
    print(f"  * Separación del peine (delta f):   {delta_f_peine:.5f} Hz")
    print(f"  * Criterio de Rayleigh (Delta f < delta f): {'CUMPLE' if df < delta_f_peine else 'NO CUMPLE (Aliasing espectral)'}")
    print(f"  * SNR recuperado sobre ruido puro:   rho = {res['SNR_rho']:.2f} sigma")

print("\n" + "=" * 85)
print("CONCLUSIÓN: Con T = 32 s hay ensanchamiento artificial de los bins.")
print("Para aislar rigurosamente el peine de 17 líneas se requieren series de T >= 256 s.")
print("=" * 85)

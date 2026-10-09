#!/usr/bin/env python3
"""
Consulta dinámica a la API de GWOSC para obtener archivos HDF5 de 4096s
en modo de ciencia limpio (O1).
"""
import os
import sys
import urllib.request
import json

print("=" * 80)
print(" LOCALIZADOR OFICIAL DE ARCHIVOS DE 4096s EN GWOSC")
print("=" * 80)

# Buscador a través de la API REST de GWOSC para la campaña O1 (4096 Hz)
# Consultamos los archivos disponibles alrededor del tiempo de referencia GPS 1126259462
t_target = 1126259462
url_api = f"https://gwosc.org/archive/links/O1/H1_GWOSC_4KHZ_R1/{t_target-4096}/{t_target+4096}/json/"

headers = {'User-Agent': 'Mozilla/5.0'}
req = urllib.request.Request(url_api, headers=headers)

try:
    print("[*] Consultando índice oficial de GWOSC...")
    with urllib.request.urlopen(req, timeout=15) as resp:
        datos = json.loads(resp.read().decode('utf-8'))
        
    archivos = datos.get('strain', [])
    print(f"[+] Total de archivos H1 localizados: {len(archivos)}")
    
    url_h1_4096 = None
    for item in archivos:
        if item.get('format') == 'hdf5' and item.get('duration') == 4096:
            url_h1_4096 = item.get('url')
            break
            
    if not url_h1_4096 and archivos:
        # Si no hay de 4096, tomar el bloque HDF5 más largo disponible
        archivos_hdf5 = [a for a in archivos if a.get('format') == 'hdf5']
        archivos_hdf5.sort(key=lambda x: x.get('duration', 0), reverse=True)
        url_h1_4096 = archivos_hdf5[0]['url']
        
    print(f"[+] URL canónica resuelta para H1:\n    {url_h1_4096}\n")
    
    # Resolver la URL correspondiente para L1
    url_l1_4096 = url_h1_4096.replace('/H1_', '/L1_').replace('H-H1', 'L-L1')
    print(f"[+] URL canónica resuelta para L1:\n    {url_l1_4096}\n")

    # Guardar las URLs en un archivo de configuración para el análisis
    with open("gwosc_urls.json", "w") as f:
        json.dump({"H1": url_h1_4096, "L1": url_l1_4096}, f, indent=2)
    print("[+] URLs guardadas con éxito en 'gwosc_urls.json'.")

except Exception as e:
    print(f"[-] Error al resolver la API: {e}")
    # Enlace canónico estático del bloque estándar que contiene GW150914 (4096 segundos)
    # Este bloque empieza en 1126257414 y dura 4096 s (GPS: 1126257414 a 1126261510)
    print("[*] Usando bloque estándar canónico de O1 (4096 s):")
    url_h1 = "https://gwosc.org/archive/data/O1/1125122048/H-H1_GWOSC_4KHZ_R1-1126256640-4096.hdf5"
    url_l1 = "https://gwosc.org/archive/data/O1/1125122048/L-L1_GWOSC_4KHZ_R1-1126256640-4096.hdf5"
    with open("gwosc_urls.json", "w") as f:
        json.dump({"H1": url_h1, "L1": url_l1}, f, indent=2)
    print("[+] URLs estáticas fijadas en 'gwosc_urls.json'.")

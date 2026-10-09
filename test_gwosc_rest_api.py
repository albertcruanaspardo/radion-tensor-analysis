#!/usr/bin/env python3
import urllib.request
import json
import gzip
import numpy as np

print("[*] Consultando API pública de GWOSC...")
# Evento GW150914 (O1)
url_api = "https://gwosc.org/eventapi/json/GWTC-1-confident/GW150914/v3/"
req = urllib.request.Request(url_api, headers={'User-Agent': 'Mozilla/5.0'})

try:
    with urllib.request.urlopen(req) as resp:
        data = json.loads(resp.read().decode('utf-8'))
        strain_files = data['events']['GW150914-v3']['strain']
        print("[+] Conexión exitosa. Archivos disponibles en GWOSC:")
        for item in strain_files:
            if item.get('format') in ['txt.gz', 'hdf5'] and item.get('sampling_rate') == 4096:
                print(f"    - Detector: {item['detector']} | Duración: {item['duration']}s | Formato: {item['format']} | URL: {item['url']}")
except Exception as e:
    print(f"[-] Error al consultar GWOSC: {e}")

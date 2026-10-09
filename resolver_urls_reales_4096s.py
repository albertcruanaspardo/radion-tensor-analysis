#!/usr/bin/env python3
"""
Consulta la API JSON oficial de GWOSC para obtener las URLs
verificadas de los archivos de 4096 segundos de H1 y L1.
"""
import urllib.request
import json

api_url = "https://gwosc.org/eventapi/json/GWTC-1-confident/GW150914/v3/"
headers = {'User-Agent': 'Mozilla/5.0'}

req = urllib.request.Request(api_url, headers=headers)
with urllib.request.urlopen(req) as resp:
    meta = json.loads(resp.read().decode('utf-8'))

strain_list = meta['events']['GW150914-v3']['strain']

urls = {}
for item in strain_list:
    det = item.get('detector')
    fmt = item.get('format')
    dur = item.get('duration')
    rate = item.get('sampling_rate')
    
    if fmt == 'hdf5' and dur == 4096 and rate == 4096:
        urls[det] = item.get('url')

print("=" * 80)
print(" ARCHIVOS HDF5 DE 4096s LOCALIZADOS EN GWOSC:")
print("=" * 80)
for det, u in urls.items():
    print(f"[{det}] {u}")

if 'H1' in urls and 'L1' in urls:
    with open("gwosc_urls.json", "w") as f:
        json.dump(urls, f, indent=2)
    print("\n[+] 'gwosc_urls.json' actualizado correctamente con enlaces verificados.")
else:
    print("\n[-] No se encontraron ambos detectores con duration=4096 y rate=4096.")
    print("Archivos disponibles con duration=4096:")
    for item in strain_list:
        if item.get('duration') == 4096:
            print(f"  Det: {item.get('detector')} | Fmt: {item.get('format')} | Rate: {item.get('sampling_rate')} | URL: {item.get('url')}")

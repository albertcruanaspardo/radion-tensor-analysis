#!/usr/bin/env bash
set -e

echo "[*] Descargando índice JSON de GW150914 con curl..."
curl -s -L --connect-timeout 10 --retry 2 \
  -H "User-Agent: Mozilla/5.0" \
  "https://gwosc.org/eventapi/json/GWTC-1-confident/GW150914/v3/" -o gw150914_meta.json

python3 - << 'PYEOF'
import json

with open("gw150914_meta.json", "r") as f:
    meta = json.load(f)

strain_list = meta['events']['GW150914-v3']['strain']

urls = {}
for item in strain_list:
    det = item.get('detector')
    fmt = item.get('format')
    dur = item.get('duration')
    rate = item.get('sampling_rate')
    
    # Prioridad: hdf5 a 4096 Hz de 4096 s o el bloque más largo disponible
    if fmt == 'hdf5' and dur == 4096 and rate == 4096:
        urls[det] = item.get('url')

print("=" * 80)
print(" ARCHIVOS HDF5 DE 4096s ENCONTRADOS:")
print("=" * 80)
for det, u in urls.items():
    print(f"[{det}] {u}")

if 'H1' in urls and 'L1' in urls:
    with open("gwosc_urls.json", "w") as f_out:
        json.dump(urls, f_out, indent=2)
    print("\n[+] gwosc_urls.json generado exitosamente.")
else:
    print("\n[-] No se encontraron ambos de 4096s a 4096Hz. Listando disponibles de 4096s:")
    for item in strain_list:
        if item.get('duration') == 4096:
            print(f"  Det: {item.get('detector')} | Fmt: {item.get('format')} | Rate: {item.get('sampling_rate')} | URL: {item.get('url')}")
PYEOF

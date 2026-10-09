#!/usr/bin/env python3
"""
Inspección de rendimientos observados vs esperados (Table 1)
de HEPData record 98796 (ATLAS Displaced Leptons 139 fb^-1, PRL 127, 051802).
"""
import urllib.request
import json

# URL oficial de la tabla 13 (Table 1: Yields) en HEPData
url = "https://www.hepdata.net/record/resource/98796.v2/t13?format=json"

headers = {
    "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko)"
}

try:
    req = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(req, timeout=10) as resp:
        table = json.loads(resp.read().decode())
        
    print("================================================================================")
    print(" DATOS REALES DE ATLAS: RENDIMIENTOS EN REGIONES DE SEÑAL (TABLE 1)")
    print(" DOI: 10.17182/hepdata.98796.v2/t13")
    print("================================================================================")
    
    # Extraer canales (SR-ee, SR-em, SR-mm)
    sr_names = [x['value'] for x in table['independent_variables'][0]['values']]
    
    # Extraer columnas dependientes (Observed, Background total, etc.)
    for dep in table['dependent_variables']:
        header = dep['header']['name']
        vals = [f"{v['value']}" for v in dep['values']]
        print(f"\nColumna: {header}")
        for sr, val in zip(sr_names, vals):
            print(f"  {sr:<10}: {val}")
            
except Exception as e:
    print(f"Acceso directo vía API no completado ({e}).")
    print("\nValores tabulados en la publicación impresa (Table 1, PRL 127, 051802):")
    print("  SR-ee:  Obs = 1 | Bkg = 1.43 +- 0.38 (stat) +- 0.35 (syst)")
    print("  SR-em:  Obs = 0 | Bkg = 0.46 +- 0.23 (stat) +- 0.17 (syst)")
    print("  SR-mm:  Obs = 1 | Bkg = 0.73 +- 0.21 (stat) +- 0.18 (syst)")
    print("  TOTAL:  Obs = 2 | Bkg = 2.62 +- 0.61")

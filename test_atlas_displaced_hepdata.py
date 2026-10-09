#!/usr/bin/env python3
import requests
import yaml
import os

record_id = "ins1841982"
api_url = f"https://www.hepdata.net/record/{record_id}?format=json"

headers = {
    "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0.0.0 Safari/537.36",
    "Accept": "application/json"
}

session = requests.Session()
response = session.get(api_url, headers=headers)

if response.status_code != 200:
    print(f"[-] Código HTTP: {response.status_code}")
    print("Cloudflare o endpoint bloqueado. Probando descarga directa de tabla YAML...")
else:
    data = response.json()
    print("================================================================================")
    print(" HEPData: REPORTE OFICIAL DEL ANÁLISIS DE ATLAS")
    print("================================================================================")
    rec = data.get("record", {})
    print(f"Título:        {rec.get('title')}")
    print(f"Colaboración:  {rec.get('collaboration')}")
    print(f"DOI/arXiv:     {rec.get('arxiv_id', 'N/A')}")
    tables = data.get("data_tables", [])
    print(f"Tablas totales digitalizadas: {len(tables)}")
    print("--------------------------------------------------------------------------------")
    print("Tablas disponibles (inspección de variables y bines):")
    for t in tables[:10]:
        print(f" -> ID {t.get('id'):2d} | Nombre: {t.get('name'):<25} | {t.get('description', '')[:65]}...")
    print("================================================================================")

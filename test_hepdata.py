#!/usr/bin/env python3
"""
Consulta directa a la API de HEPData (CERN) para inspeccionar datos
reales de búsquedas de vértices desplazados (Displaced Vertices / LLPs).
"""
import urllib.request
import json

record_id = "ins1841982"
url = f"https://www.hepdata.net/record/{record_id}?format=json"

print(f"Consultando HEPData para el registro {record_id}...")

try:
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req, timeout=15) as response:
        data = json.loads(response.read().decode('utf-8'))
        
        print("\n--- DATOS DE LA PUBLICACIÓN OFICIAL ---")
        rec = data.get('record', {})
        print(f"Título: {rec.get('title', 'N/A')}")
        print(f"Colaboración: {rec.get('collaboration', 'N/A')}")
        tables = data.get('data_tables', [])
        print(f"Número de tablas disponibles: {len(tables)}")
        
        print("\nPrimeras tablas disponibles:")
        for t in tables[:8]:
            print(f" - ID {t.get('id')}: {t.get('name')} | {t.get('description', '')[:70]}...")
            
except Exception as e:
    print(f"No se pudo consultar el registro o no hay conexión: {e}")

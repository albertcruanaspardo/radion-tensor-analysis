#!/usr/bin/env python3
import urllib.request
import json

url = "https://inspirehep.net/api/literature?q=arxiv:2011.07812"

req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
with urllib.request.urlopen(req) as resp:
    data = json.loads(resp.read().decode())

hits = data.get("hits", {}).get("hits", [])
if not hits:
    print("[-] No se encontró el registro en INSPIRE.")
    exit(1)

metadata = hits[0].get("metadata", {})
print("================================================================================")
print(" PUBLICACIÓN OFICIAL DE ATLAS")
print("================================================================================")
print("Título:", metadata.get("titles", [{}])[0].get("title"))
print("Publicación:", metadata.get("publication_info", [{}])[0].get("journal_title"),
      metadata.get("publication_info", [{}])[0].get("journal_volume"),
      f"({metadata.get('publication_info', [{}])[0].get('year')})")
print("DOI:", metadata.get("dois", [{}])[0].get("value"))

print("\nEnlaces a datos suplementarios y repositorios:")
for url_entry in metadata.get("urls", []):
    print(" ->", url_entry.get("value"))

for ext in metadata.get("external_system_identifiers", []):
    print(f" -> {ext.get('schema')}: {ext.get('value')}")
print("================================================================================")

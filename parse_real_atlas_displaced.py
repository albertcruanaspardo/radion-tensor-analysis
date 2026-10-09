#!/usr/bin/env python3
import os
import tarfile
import yaml

archive = "ins1831504.tar.gz"

if not os.path.exists(archive) or not tarfile.is_tarfile(archive):
    print(f"[-] No se encontró un tarball válido: {archive}")
    print("Descárgalo desde: https://www.hepdata.net/record/ins1831504")
    exit(1)

out_dir = "atlas_ins1831504_raw"
os.makedirs(out_dir, exist_ok=True)
with tarfile.open(archive, "r:gz") as tar:
    tar.extractall(path=out_dir)

sub_file = os.path.join(out_dir, "submission.yaml")
with open(sub_file, "r") as f:
    docs = list(yaml.safe_load_all(f))

print("================================================================================")
print(" TABLAS DE DATOS REALES DE ATLAS (PRL 127, 051802)")
print("================================================================================")

for doc in docs:
    if not doc or "name" not in doc:
        continue
    name = doc.get("name")
    data_file = doc.get("data_file")
    desc = doc.get("description", "")
    print(f"\n[Tabla] {name} -> Archivo: {data_file}")
    print(f"  Detalle: {desc[:100]}...")
    
    if data_file:
        table_path = os.path.join(out_dir, data_file)
        if os.path.exists(table_path):
            with open(table_path, "r") as tf:
                tdata = yaml.safe_load(tf)
                indep = [v.get("header", {}).get("name") for v in tdata.get("independent_variables", [])]
                dep = [v.get("header", {}).get("name") for v in tdata.get("dependent_variables", [])]
                print(f"  X (Bines): {indep}")
                print(f"  Y (Mediciones): {dep}")

print("================================================================================")

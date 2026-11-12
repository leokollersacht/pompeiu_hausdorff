#!/usr/bin/env python3
"""
Automated benchmark and evidence generation for Checkpoint 1.
Project: Pompeiu-Hausdorff distance (Sacht & Jacobson, 2024)
"""

import os
import sys
import subprocess
import time
from datetime import datetime

def run_baseline():
    repo_dir = os.path.dirname(os.path.abspath(__file__))
    build_dir = os.path.join(repo_dir, "build")
    executable = os.path.join(build_dir, "pompeiu_hausdorff")
    mesh_a = os.path.join(repo_dir, "meshes", "107100.obj")
    mesh_b = os.path.join(repo_dir, "meshes", "107100_sf.obj")
    output_file = os.path.join(repo_dir, "baseline_output.txt")
    
    tol = "1e-4"
    max_factor = "1000000"
    normalize = "1"

    if not os.path.isfile(executable):
        print(f"Erro: Binário não encontrado em {executable}.")
        print("Certifique-se de compilar o projeto antes de executar o script.")
        sys.exit(1)

    if not os.path.isfile(mesh_a) or not os.path.isfile(mesh_b):
        print(f"Erro: Malhas de teste não encontradas em {mesh_a} ou {mesh_b}.")
        sys.exit(1)

    cmd = [executable, mesh_a, mesh_b, tol, max_factor, normalize]
    print(f"Executando baseline:")
    print(" ".join(cmd))
    print("-" * 50)

    start_time = time.time()
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    elapsed = time.time() - start_time

    timestamp_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    report = []
    report.append("=" * 60)
    report.append("CHECKPOINT 1 - EVIDÊNCIA DE EXECUÇÃO DE BASELINE")
    report.append("Artigo: Cascading upper bounds for triangle soup Pompeiu-Hausdorff distance")
    report.append(f"Data e Hora: {timestamp_str}")
    report.append(f"Comando executado: {' '.join(cmd)}")
    report.append(f"Tempo total de parede (wall-clock): {elapsed:.4f}s")
    report.append(f"Código de retorno do processo: {result.returncode}")
    report.append("=" * 60)
    report.append("\n--- [SAÍDA PADRÃO (STDOUT)] ---")
    report.append(result.stdout.strip())
    if result.stderr.strip():
        report.append("\n--- [ERRO PADRÃO (STDERR)] ---")
        report.append(result.stderr.strip())
    report.append("\n" + "=" * 60)

    content = "\n".join(report) + "\n"

    with open(output_file, "w", encoding="utf-8") as f:
        f.write(content)

    print(content)
    print(f"[OK] Saída salva com sucesso em: {output_file}")

if __name__ == "__main__":
    run_baseline()

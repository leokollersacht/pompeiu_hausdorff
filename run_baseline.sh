#!/usr/bin/env bash
# Script Bash para execução e registro do Checkpoint 1 (Baseline)
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="$SCRIPT_DIR/build/pompeiu_hausdorff"
MESH_A="$SCRIPT_DIR/meshes/107100.obj"
MESH_B="$SCRIPT_DIR/meshes/107100_sf.obj"
OUTPUT_FILE="$SCRIPT_DIR/baseline_output.txt"
TOL="1e-4"
MAX_FACTOR="1000000"
NORMALIZE="1"

if [ ! -f "$BIN" ]; then
    echo "Erro: Executável não encontrado em $BIN. Compile o projeto primeiro."
    exit 1
fi

echo "Executando baseline original..."
{
    echo "============================================================"
    echo "CHECKPOINT 1 - EVIDÊNCIA DE EXECUÇÃO DE BASELINE"
    echo "Data: $(date)"
    echo "Comando: $BIN $MESH_A $MESH_B $TOL $MAX_FACTOR $NORMALIZE"
    echo "============================================================"
    echo ""
    # Executa o binário (o código original retorna 1 em caso de sucesso, por isso '|| true')
    "$BIN" "$MESH_A" "$MESH_B" "$TOL" "$MAX_FACTOR" "$NORMALIZE" || true
    echo ""
    echo "============================================================"
} | tee "$OUTPUT_FILE"

echo "[OK] Evidência gravada com sucesso em: $OUTPUT_FILE"

#!/usr/bin/env bash
# Paso 3. Elementos moviles y grupos biosinteticos.
#
# NOTA IMPORTANTE SOBRE geNomad: clasifica plasmidos y virus. NO delimita
# elementos conjugativos integrativos, y su conjunto de marcadores anota mal
# fuera de ese dominio (en esta corrida llamo "recombinasa" a metG). Los
# limites de los dos ICE del manuscrito se establecieron a partir de la
# anotacion de Bakta, gen por gen, no de la salida de geNomad.
set -euo pipefail
source "$(dirname "$0")/config.sh"
cd "$DATA_DIR"

run "$ENV_MGE" genomad end-to-end \
    flye/S1/assembly.fasta genomad_S1 "$GENOMAD_DB" \
    --threads "$THREADS"

# antiSMASH exige un directorio de salida vacio
rm -rf antismash_S1
run "$ENV_AS" antismash \
    --output-dir antismash_S1 \
    --genefinding-tool prodigal \
    --cb-knownclusters --cb-subclusters \
    --cpus "$THREADS" \
    flye/S1/assembly.fasta

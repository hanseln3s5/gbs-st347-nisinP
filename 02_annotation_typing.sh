#!/usr/bin/env bash
# Paso 2. Anotacion, tipificacion y determinantes de resistencia del aislado.
set -euo pipefail
source "$(dirname "$0")/config.sh"
cd "$DATA_DIR"

# Anotacion estructural y funcional
run "$ENV_ANNOT" bakta \
    --db "$BAKTA_DB" \
    --output bakta_S1 \
    --prefix S1 \
    --threads "$THREADS" \
    flye/S1/assembly.fasta

# MLST (esquema de S. agalactiae: adhP pheS atr glnA sdhA glcK tkt)
run "$ENV_TYPING" mlst flye/S1/assembly.fasta flye/S2/assembly.fasta \
    > results_mlst.tsv

# Serotipo capsular. GBS-SBG necesita blastn en el PATH: por eso va en
# gbs-typing y no en gbs-mge, donde blastn no existe.
run "$ENV_TYPING" bash -c \
    'python GBS-SBG/GBS-SBG.py -i flye/S1/assembly.fasta' \
    > results_serotipo_S1.txt

# Determinantes de resistencia. Se corre con las DOS versiones de base de
# datos porque el manuscrito reporta concordancia entre ellas.
run "$ENV_TYPING" amrfinder -n flye/S1/assembly.fasta --plus \
    --organism Streptococcus_agalactiae > results_amrfinder_db2024.tsv
run "$ENV_ANNOT"  amrfinder -n flye/S1/assembly.fasta --plus \
    --organism Streptococcus_agalactiae > results_amrfinder_db2026.tsv

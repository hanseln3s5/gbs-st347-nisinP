#!/usr/bin/env bash
# Paso 5. Comparacion de las dos librerias.
set -euo pipefail
source "$(dirname "$0")/config.sh"
cd "$DATA_DIR"

# Alineamiento de los dos ensamblajes
run "$ENV_CMP" dnadiff flye/S1/assembly.fasta flye/S2/assembly.fasta -p s1s2

# Fracciones alelicas en las posiciones discordantes y prueba de homogeneidad.
# El script localiza show-snps (gbs-cmp) y samtools (gbs-asm) por si mismo.
run "$ENV_ASM" python3 scripts/fracciones.py \
    flye/S1/assembly.fasta S1.bam S2_vs_S1.bam \
    | tee results_fracciones_alelicas.txt

# Control de particion: ¿puede una particion aleatoria de un solo pool
# producir la diferencia de N50 observada entre S1 y S2?
bash scripts/particion_control.sh

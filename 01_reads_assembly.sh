#!/usr/bin/env bash
# Paso 1. Estadisticas de lecturas, filtrado y ensamblaje de S1 y S2.
set -euo pipefail
source "$(dirname "$0")/config.sh"
cd "$DATA_DIR"

mkdir -p filt flye

for S in S1 S2; do
  # estadisticas de las lecturas crudas
  run "$ENV_ASM" nanoq -i "${S}.fastq" -s > "results_${S}_raw_stats.txt"

  # filtrado
  run "$ENV_ASM" filtlong \
      --min_length 1000 \
      --keep_percent 90 \
      --target_bases 420000000 \
      "${S}.fastq" | gzip > "filt/${S}.filt.fastq.gz"

  run "$ENV_ASM" nanoq -i "filt/${S}.filt.fastq.gz" -s > "results_${S}_filt_stats.txt"

  # ensamblaje
  run "$ENV_ASM" flye \
      --nano-hq "filt/${S}.filt.fastq.gz" \
      --out-dir "flye/${S}" \
      --threads "$THREADS"
done

# mapeo de cada librería contra el ensamblaje de S1 (coordenadas comunes)
run "$ENV_ASM" minimap2 -ax map-ont -t "$THREADS" \
    flye/S1/assembly.fasta filt/S1.filt.fastq.gz \
  | run "$ENV_ASM" samtools sort -@ 8 -o S1.bam -
run "$ENV_ASM" samtools index S1.bam

run "$ENV_ASM" minimap2 -ax map-ont -t "$THREADS" \
    flye/S1/assembly.fasta filt/S2.filt.fastq.gz \
  | run "$ENV_ASM" samtools sort -@ 8 -o S2_vs_S1.bam -
run "$ENV_ASM" samtools index S2_vs_S1.bam

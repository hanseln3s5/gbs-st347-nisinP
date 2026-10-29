#!/usr/bin/env bash
# Paso 4. Tamizaje de los 83 genomas dominicanos publicos.
# Los genomas provienen de BioProject PRJNA945321 (Laycock et al. 2023);
# sus accesiones estan en data/accessions.txt.
set -euo pipefail
source "$(dirname "$0")/config.sh"
cd "$DATA_DIR"

# tblastn de las proteinas consulta contra cada genoma
run "$ENV_TYPING" bash -c '
  for f in dr_genomes/*.fna; do
    acc=$(basename "$f" .fna)
    tblastn -query clusterP_prot.faa -subject "$f" \
            -outfmt "6 qseqid sseqid pident length evalue bitscore" \
            -evalue 1e-10 \
      | awk -v a="$acc" "{print \$1\"\t\"a\"\t\"\$3\"\t\"\$5}"
  done
' > operonP_dr.tsv

run "$ENV_TYPING" bash -c '
  for f in dr_genomes/*.fna; do
    acc=$(basename "$f" .fna)
    tblastn -query amr_query.faa -subject "$f" \
            -outfmt "6 qseqid sseqid pident length evalue bitscore" \
            -evalue 1e-10 \
      | awk -v a="$acc" "{print a\"\t\"\$1\"\t\"\$3\"\t\"\$5}"
  done
' > amr_dr.tsv

# Resistoma completo con AMRFinderPlus
run "$ENV_ANNOT" bash -c '
  for f in dr_genomes/*.fna; do
    acc=$(basename "$f" .fna)
    amrfinder -n "$f" --plus --organism Streptococcus_agalactiae \
      | awk -v a="$acc" "NR>1{print a\"\t\"\$6}"
  done
' > resistoma_dr.tsv

# ST y serotipo de cada genoma
run "$ENV_TYPING" bash -c '
  for f in dr_genomes/*.fna; do
    acc=$(basename "$f" .fna)
    st=$(mlst "$f" | cut -f3)
    sero=$(python GBS-SBG/GBS-SBG.py -i "$f" | tail -1 | cut -f2)
    printf "%s\t%s\t%s\n" "$acc" "$st" "$sero"
  done
' | sed 's/III-[1-4]/III/' > st_serotipo_dr.tsv

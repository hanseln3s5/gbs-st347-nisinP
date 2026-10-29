#!/usr/bin/env bash
# Variables compartidas por todos los pasos. Editar SOLO este archivo.
# Se usa `set -u`, asi que una variable sin definir detiene el script
# en vez de producir una ruta vacia silenciosa.

# Directorio de trabajo con los datos del proyecto
DATA_DIR="${DATA_DIR:-/mnt/das2/project2/data}"

# Hilos
THREADS="${THREADS:-16}"

# Bases de datos locales (rutas en la maquina donde se corrio)
BAKTA_DB="${BAKTA_DB:?definir la ruta a la base de datos de Bakta}"
GENOMAD_DB="${GENOMAD_DB:?definir la ruta a la base de datos de geNomad}"

# Entornos de micromamba usados en cada paso
ENV_ASM=gbs-asm        # filtlong, flye, minimap2, nanoq, samtools, seqkit, scipy
ENV_TYPING=gbs-typing  # blast, mlst, amrfinderplus (DB 2024-07-22.1)
ENV_ANNOT=gbs-annot    # bakta, amrfinderplus (DB 2026-08-07.1)
ENV_MGE=gbs-mge        # genomad
ENV_AS=antismash       # antismash 8
ENV_CMP=gbs-cmp        # mummer4 (dnadiff, show-snps), seqkit

run() { micromamba run -n "$1" "${@:2}"; }

# Datos

Este directorio guarda las **tablas derivadas** del analisis, no los datos crudos.
Las lecturas, los ensamblajes, los BAM y las bases de datos quedan fuera del
repositorio (ver `.gitignore`): son grandes y, en el caso de las lecturas, van a
la SRA.

## Lo que va aqui

Copiar desde el directorio de trabajo una vez corrido el pipeline:

| Archivo | Lo produce | Contenido |
|---|---|---|
| `st_serotipo_dr.tsv` | paso 4 | accesion, ST y serotipo de cada genoma |
| `operonP_dr.tsv` | paso 4 | hits de tblastn de las proteinas del operon |
| `amr_dr.tsv` | paso 4 | hits de tblastn de las proteinas de resistencia |
| `resistoma_dr.tsv` | paso 4 | genes de resistencia por genoma (AMRFinderPlus) |
| `Tabla_S1_genomas.tsv` | `generar_suplementarios.sh` | un genoma por fila, la entrada de `estadisticas.py` |
| `accessions.txt` | ver abajo | las accesiones de los 83 genomas de comparacion |
| `results_fracciones_alelicas.txt` | paso 5 | conteos alelicos en las posiciones discordantes |

## accessions.txt

Los 83 genomas provienen de BioProject **PRJNA945321** (Laycock et al. 2023,
PLOS Glob Public Health 3(9):e0002281). El estudio original secuencio 93
aislados; 83 ensamblajes estaban disponibles para descarga, y esa diferencia es
una limitacion de este reanalisis, no del trabajo original.

Para generar la lista desde el directorio de genomas:

    ls dr_genomes/*.fna | xargs -n1 basename | sed 's/\.fna$//' > data/accessions.txt

## Lo que NO va aqui

Nada que identifique a una persona. El aislado secuenciado procede de material
residual anonimizado y no se transfirio informacion clinica ni demografica.

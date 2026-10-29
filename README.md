# A closed genome of *Streptococcus agalactiae* ST347/II: a composite integrative element coupling nisin P biosynthesis to *erm*(B) and *tet*(O)

Code and derived data for the manuscript:

> García-Encarnación L, De La Cruz-Caminero DS, López-Chevalier HB.
> *A closed genome of Streptococcus agalactiae ST347/II reveals a composite
> integrative element coupling nisin P biosynthesis to erm(B) and tet(O) in
> Dominican group B Streptococcus.* Preprint, 2026.
> doi: **[PENDING — bioRxiv DOI]**

Correspondence: Hansel B. López-Chevalier (hlopez@ipl.edu.do), Instituto
Especializado de Estudios Superiores Loyola, San Cristóbal, Dominican Republic.

---

## What this repository contains

Everything needed to reproduce the analysis from the raw reads and the public
comparison genomes: the pipeline as it was actually run, the statistics, the
figure scripts, and the derived tables the manuscript reports.

It does **not** contain raw reads, assemblies or databases. Those go to the SRA
and GenBank (see *Data availability*), and the 83 comparison genomes are public
and are listed by accession rather than redistributed.

```
pipeline/     the five steps, in order, with the exact parameters used
scripts/      analysis and figure scripts called by the pipeline
data/         derived tables (see data/README.md for what belongs here)
environment/  conda/micromamba environment specifications
figures/      the three manuscript figures as rendered
results/      output directory, empty in the repository
```

## Quick start

```bash
git clone https://github.com/<user>/gbs-st347-nisinP.git
cd gbs-st347-nisinP

# 1. recreate the environments (see environment/README.md)
micromamba env create -f environment/gbs-asm.yml

# 2. point the pipeline at your data and databases
$EDITOR pipeline/config.sh

# 3. run the steps in order
bash pipeline/01_reads_assembly.sh
bash pipeline/02_annotation_typing.sh
bash pipeline/03_mobile_elements.sh
bash pipeline/04_comparative_screen.sh
bash pipeline/05_replicates.sh

# 4. build the supplementary tables and the statistics
bash scripts/generar_suplementarios.sh
python3 scripts/estadisticas.py suplementarios/Tabla_S1_genomas.tsv
python3 scripts/figuras.py
```

`pipeline/config.sh` is the only file you should need to edit. It uses `set -u`
semantics, so an undefined database path stops the run instead of silently
producing an empty one.

## The analysis, step by step

**1. Reads and assembly.** Read statistics with `nanoq`; filtering with
`filtlong --min_length 1000 --keep_percent 90 --target_bases 420000000`
(≈203× coverage of the 2.07 Mb chromosome); assembly with `flye --nano-hq`, the
mode intended for R10.4.1 chemistry with super-accuracy basecalling. Both
libraries are then mapped to the S1 assembly so that every later measurement is
in one coordinate system.

**2. Annotation and typing.** `bakta` for annotation; `mlst` for the
*S. agalactiae* scheme (*adhP*, *pheS*, *atr*, *glnA*, *sdhA*, *glcK*, *tkt*);
`GBS-SBG` for capsular serotype, with serotype III subtype calls (III-1 to
III-4) collapsed to III; `AMRFinderPlus` under two database versions, because
the manuscript reports concordance between them.

**3. Mobile elements and biosynthetic clusters.** `geNomad` and `antiSMASH`.

> **A caution recorded here because it cost us time.** geNomad classifies
> plasmids and viruses. It does **not** delimit integrative conjugative
> elements, and its marker set misannotates outside that domain — in this run
> it called *metG* a recombinase, which led us to a wrong conclusion about one
> element's mobility until Bakta's annotation corrected it. The boundaries of
> both ICEs in the manuscript were established gene by gene from Bakta, not
> from geNomad's output.

**4. Comparative screen.** `tblastn` of the cluster and resistance query
proteins against each of the 83 Dominican genomes, plus `AMRFinderPlus` for the
full resistome, and `mlst` and `GBS-SBG` for lineage and serotype.

> Operon presence is defined by **seven** cluster-specific genes. Three query
> proteins — *nisT*, *nisR* and the ABC ATPase — are excluded from the
> criterion because they cross-hit native streptococcal transporters and
> regulators in 81 to 83 of the 83 genomes and therefore carry no information.
> A different gene set would give a different count; the criterion is stated so
> that it can be reproduced or disputed.

**5. Replicate comparison.** `dnadiff` to align the two assemblies, then
`scripts/fracciones.py` to tabulate allele counts at the discordant positions
and test whether a single mixing ratio accounts for them, and
`scripts/particion_control.sh` to ask whether a random partition of one read
pool reproduces the observed difference in read N50.

## Scripts

| Script | What it does |
|---|---|
| `scripts/fracciones.py` | Allele fractions at the positions where the two assemblies differ, with exact binomial confidence intervals and a chi-square test of homogeneity across positions. Locates `show-snps` and `samtools` across micromamba environments by itself. |
| `scripts/ligamiento.py` | Attempts read-level phasing of the discordant positions. **Kept for the record; it cannot work on this dataset** — the three substitutions lie 260 kb, 596 kb and 856 kb apart and no read spans a pair. It is included because the manuscript states that phasing was attempted and why it fails. |
| `scripts/estadisticas.py` | Fisher's exact tests and odds ratios (Haldane-Anscombe corrected where a cell is zero) for the association between the operon and the resistance determinants, plus the lineage breakdowns. |
| `scripts/particion_control.sh` | The empirical partition control. |
| `scripts/generar_suplementarios.sh` | Builds the supplementary tables from the pipeline output. |
| `scripts/figuras.py` | Figures 2 and 3. |

## Key results these scripts reproduce

- The operon is present in **12 of 83** comparison genomes (13 of 84, 15.5%,
  including the sequenced isolate), with an invariant 55-residue precursor.
- *tet*(O) in **12/12** carriers against **0/71** others, Fisher
  p = 1.03 × 10⁻¹⁴; *erm*(B) in **12/12** against **3/71**, p = 4.70 × 10⁻¹².
- The nisin resistance protein is present in **72/83**; the eleven absences are
  the entire ST1 lineage, an exact lineage split and a testable prediction of
  susceptibility.
- All **15** ST1010/IV genomes carry *aac*(6′)-Ie/*aph*(2″)-Ia and no other
  resistance determinant.
- Random partitions of a single read pool give read N50 of 7,463 / 7,461 /
  7,457 bp — a 6 bp spread against the 392 bp observed between S1 and S2.
- Minor-allele fractions at the three discordant substitutions are 29.4%,
  46.6% and 55.8% pooled, heterogeneous at χ² = 59.25, 2 df, p = 1.4 × 10⁻¹³,
  which excludes a fixed-ratio two-clone mixture.

## Data availability

| What | Where |
|---|---|
| Raw reads (two runs, one BioSample) | SRA — **[PENDING accession]** |
| Closed assembly | GenBank — **[PENDING accession]** |
| 83 comparison genomes | BioProject **PRJNA945321** (Laycock et al. 2023); accessions in `data/accessions.txt` |
| Derived tables | `data/` in this repository |

The 83 comparison genomes are **maternal colonisation isolates** from a
published Dominican survey of near-term pregnant women — not clinical or
invasive disease isolates. The sequenced isolate comes from a different
sampling frame: routine diagnostic material at a reference laboratory.

## Licence

Code (`scripts/`, `pipeline/`) under the MIT licence; see `LICENSE`.
Derived data (`data/`) and figures under CC BY 4.0; see `LICENSE-DATA`.

## Citing

See `CITATION.cff`. 

## Reproducibility notes

- Environment specifications in `environment/` must be **exported from the
  machine that ran the analysis**, not written by hand, or the versions will
  drift from what was used.
- Two AMRFinderPlus database versions were used deliberately (2024-07-22.1 and
  2026-08-07.1). Results concord; re-running with only one is fine but will not
  reproduce the concordance statement.
- `mlst` ST numbers are registry indices assigned in order of submission to
  PubMLST. They carry no biological meaning on their own; clonal complex
  assignment (sharing ≥5/7 alleles with a founder) is what groups lineages.

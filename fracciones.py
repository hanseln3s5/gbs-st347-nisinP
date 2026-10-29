#!/usr/bin/env python3
"""
Fracciones alelicas en las posiciones que difieren entre S1 y S2.

Prueba la hipotesis de dos clones co-residentes: una mezcla de dos clones
tiene UNA proporcion, asi que todos los sitios que los distingan deben
mostrar la misma fraccion. Si las fracciones son heterogeneas, no hay
mezcla a proporcion fija y el patron es error sistematico del basecaller.

Uso:
    python3 fracciones.py REF.fasta BAM [BAM2 ...]

Ejemplo:
    python3 fracciones.py flye/S1/assembly.fasta S1.bam S2_vs_S1.bam

No toma posiciones como argumento: las detecta del .delta si existe,
o usa las tres ya conocidas.
"""

import collections
import glob
import os
import shutil
import subprocess
import sys

POSICIONES_CONOCIDAS = [1198656, 1794246, 2054407]


def tool(name):
    """Encuentra un binario en PATH o en cualquier entorno de micromamba."""
    p = shutil.which(name)
    if p:
        return p
    pats = [os.path.expanduser('~/micromamba/envs/*/bin/' + name),
            os.path.expanduser('~/.micromamba/envs/*/bin/' + name),
            '/opt/*/envs/*/bin/' + name]
    for pat in pats:
        hits = sorted(glob.glob(pat))
        if hits:
            return hits[0]
    sys.exit('No encuentro %s. MUMmer esta en gbs-cmp, samtools en gbs-asm.'
             % name)


SAMTOOLS = tool('samtools')


def parse_bases(s):
    """Parsea la columna de bases de mpileup, saltando marcadores."""
    out, i, n = [], 0, len(s)
    while i < n:
        c = s[i]
        if c == '^':
            i += 2
        elif c == '$':
            i += 1
        elif c in '+-':
            j, num = i + 1, ''
            while j < n and s[j].isdigit():
                num += s[j]
                j += 1
            i = j + (int(num) if num else 0)
        elif c in '.,':
            out.append('REF')
            i += 1
        elif c.upper() in 'ACGTN':
            out.append(c.upper())
            i += 1
        elif c == '*':
            out.append('DEL')
            i += 1
        else:
            i += 1
    return out


def detectar_posiciones():
    """Lee las sustituciones del .delta si hay uno a mano."""
    deltas = sorted(glob.glob('*.delta'))
    if not deltas:
        return None, None
    showsnps = shutil.which('show-snps')
    if not showsnps:
        hits = sorted(glob.glob(
            os.path.expanduser('~/micromamba/envs/*/bin/show-snps')))
        if not hits:
            return None, None
        showsnps = hits[0]
    out = subprocess.run([showsnps, '-ClrTH', deltas[0]],
                         capture_output=True, text=True)
    if out.returncode != 0:
        return None, None
    pos, contig = [], None
    for line in out.stdout.splitlines():
        f = line.split('\t')
        if len(f) < 12:
            continue
        p, rb, qb = f[0], f[1], f[2]
        if rb in 'ACGT' and qb in 'ACGT':      # sustitucion, no indel
            pos.append(int(p))
            contig = f[10]
    return (sorted(set(pos)), contig) if pos else (None, None)


def fracciones(ref, bam, contig, posiciones):
    """Cuenta bases en cada posicion de un BAM."""
    region = ','.join('%s:%d-%d' % (contig, p, p) for p in posiciones)
    cmd = [SAMTOOLS, 'mpileup', '-f', ref, '-aa', '-Q', '0', '-q', '0',
           '-d', '10000', '-r', region.split(',')[0], bam]
    # mpileup -r toma una sola region; iteramos
    resultado = {}
    for p in posiciones:
        cmd = [SAMTOOLS, 'mpileup', '-f', ref, '-Q', '0', '-q', '0',
               '-d', '10000', '-r', '%s:%d-%d' % (contig, p, p), bam]
        out = subprocess.run(cmd, capture_output=True, text=True)
        for line in out.stdout.splitlines():
            f = line.split('\t')
            if len(f) < 5 or int(f[1]) != p:
                continue
            refbase = f[2].upper()
            bases = parse_bases(f[4])
            cnt = collections.Counter(
                refbase if b == 'REF' else b for b in bases)
            resultado[p] = (refbase, cnt)
    return resultado


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    ref, bams = sys.argv[1], sys.argv[2:]
    for f in [ref] + bams:
        if not os.path.exists(f):
            sys.exit('No existe: %s' % f)

    posiciones, contig = detectar_posiciones()
    if not posiciones:
        posiciones = POSICIONES_CONOCIDAS
        out = subprocess.run([SAMTOOLS, 'idxstats', bams[0]],
                             capture_output=True, text=True)
        contig = out.stdout.split('\t')[0]
        print('(sin .delta legible: uso las posiciones conocidas)')
    print('contig: %s   posiciones: %s\n'
          % (contig, ', '.join(str(p) for p in posiciones)))

    tabla = {}
    for bam in bams:
        res = fracciones(ref, bam, contig, posiciones)
        print('=== %s ===' % bam)
        filas = []
        for p in posiciones:
            if p not in res:
                print('  %d  sin cobertura' % p)
                continue
            refbase, cnt = res[p]
            prof = sum(cnt.values())
            # los dos alelos mas frecuentes
            top = cnt.most_common()
            principal = refbase if refbase in cnt else top[0][0]
            alt = next((b for b, _ in top if b != principal), None)
            n_ref, n_alt = cnt.get(principal, 0), cnt.get(alt, 0) if alt else 0
            n = n_ref + n_alt
            detalle = '  '.join('%s:%d (%.0f%%)' % (b, c, 100.0 * c / prof)
                                for b, c in top)
            print('  %d  ref=%s  prof=%d  %s' % (p, refbase, prof, detalle))
            if n:
                print('      fraccion alternativa (%s) = %.1f%%  [%d/%d]'
                      % (alt, 100.0 * n_alt / n, n_alt, n))
                filas.append((p, n_ref, n_alt))
        tabla[bam] = filas
        print()

    # prueba de homogeneidad por BAM
    try:
        from scipy.stats import binomtest, chi2_contingency
    except ImportError:
        print('scipy no disponible: omito la prueba de homogeneidad.')
        return

    for bam, filas in tabla.items():
        if len(filas) < 2:
            continue
        print('--- homogeneidad de fracciones, %s ---' % bam)
        for p, r, a in filas:
            bt = binomtest(a, r + a, 0.5)
            ci = bt.proportion_ci(confidence_level=0.95)
            print('  %d  %.1f%%  IC95%% %.1f-%.1f%%  p(vs 50%%)=%.3g'
                  % (p, 100.0 * a / (r + a), 100 * ci.low, 100 * ci.high,
                     bt.pvalue))
        chi2, pval, dof, _ = chi2_contingency(
            [[r for _, r, _ in filas], [a for _, _, a in filas]])
        print('  chi2 = %.2f   gl = %d   p = %.3g' % (chi2, dof, pval))
        if pval < 0.05:
            print('  -> fracciones HETEROGENEAS: no hay una proporcion de')
            print('     mezcla unica. Incompatible con dos clones a')
            print('     proporcion fija; compatible con error sistematico.')
        else:
            print('  -> fracciones compatibles con una proporcion comun:')
            print('     NO se puede descartar una mezcla de dos clones.')
        print()


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Decide si las posiciones discordantes entre S1 y S2 estan ligadas.

Dos clones => los alelos co-segregan: las lecturas llevan una de dos
combinaciones. Error de basecalling => las posiciones son independientes y
aparecen las cuatro combinaciones en proporciones parecidas.

Uso:  python3 ligamiento.py [BAM] [REFERENCIA] [DELTA]
"""
import collections, glob, os, shutil, subprocess, sys


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
    sys.exit('No encuentro %s. Instalalo o activa el entorno que lo tiene '
             '(MUMmer esta en gbs-cmp, samtools en gbs-asm).' % name)


SHOWSNPS = tool('show-snps')
SAMTOOLS = tool('samtools')

BAM = sys.argv[1] if len(sys.argv) > 1 else 'S1.bam'
REF = sys.argv[2] if len(sys.argv) > 2 else 'flye/S1/assembly.fasta'
DELTA = sys.argv[3] if len(sys.argv) > 3 else None


def find_delta():
    if DELTA and os.path.exists(DELTA):
        return DELTA
    for pat in ('*.delta', '*/*.delta'):
        hits = [f for f in glob.glob(pat) if not f.endswith('.1delta')]
        if hits:
            return sorted(hits, key=os.path.getmtime)[-1]
    return None


def substitutions(delta):
    """Posiciones de sustitucion (no indels) en el alineamiento de replicas."""
    out = subprocess.run([SHOWSNPS, '-ClrTH', delta],
                         capture_output=True, text=True)
    if out.returncode != 0:
        sys.exit('show-snps fallo:\n' + out.stderr.strip())
    pos = []
    for line in out.stdout.splitlines():
        f = line.split('\t')
        if len(f) < 4:
            continue
        p, a, b = f[0], f[1], f[2]
        if a != '.' and b != '.':
            pos.append((int(p), a, b))
    return pos


def parse_bases(s):
    """Una entrada por lectura, en el mismo orden que la lista de QNAME."""
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
                num += s[j]; j += 1
            i = j + (int(num) if num else 0)
        elif c in '.,':
            out.append('REF'); i += 1
        elif c.upper() in 'ACGTN':
            out.append(c.upper()); i += 1
        elif c == '*':
            out.append('DEL'); i += 1
        else:
            i += 1
    return out


def column(contig, p):
    r = subprocess.run(
        [SAMTOOLS, 'mpileup', '-f', REF, '-r', '%s:%d-%d' % (contig, p, p),
         '-Q', '0', '-q', '0', '--output-QNAME', BAM],
        capture_output=True, text=True)
    line = r.stdout.strip()
    if not line:
        return None, {}
    f = line.split('\t')
    if len(f) < 7:
        return None, {}
    ref = f[2].upper()
    calls = parse_bases(f[4])
    names = f[6].strip().split(',')
    if len(calls) != len(names):
        print('  aviso: %d bases contra %d nombres en la posicion %d' %
              (len(calls), len(names), p))
    d = {}
    for name, c in zip(names, calls):
        d[name] = ref if c == 'REF' else c
    return ref, d


def contig_name():
    r = subprocess.run([SAMTOOLS, 'idxstats', BAM], capture_output=True, text=True)
    return r.stdout.split('\t')[0]


def main():
    delta = find_delta()
    if not delta:
        sys.exit('No encuentro el archivo .delta. Corre primero:\n'
                 '  dnadiff flye/S1/assembly.fasta flye/S2/assembly.fasta -p rep')
    print('delta:', delta, '  bam:', BAM)
    print('show-snps:', SHOWSNPS)
    print('samtools: ', SAMTOOLS)
    subs = substitutions(delta)
    if len(subs) < 2:
        sys.exit('Solo %d sustitucion(es); hacen falta dos o mas para medir '
                 'ligamiento.' % len(subs))
    print('sustituciones:', ', '.join('%d (%s/%s)' % s for s in subs))
    ctg = contig_name()
    print('contig:', ctg)
    print()

    cols = {}
    for p, _, _ in subs:
        ref, d = column(ctg, p)
        cols[p] = d
        if d:
            c = collections.Counter(d.values())
            tot = sum(c.values())
            print('posicion %d  ref=%s  profundidad=%d  %s' % (
                p, ref, tot,
                '  '.join('%s:%d (%.0f%%)' % (b, n, 100.0 * n / tot)
                          for b, n in c.most_common())))
    print()

    algo = False
    for i in range(len(subs)):
        for j in range(i + 1, len(subs)):
            pa, pb = subs[i][0], subs[j][0]
            shared = set(cols[pa]) & set(cols[pb])
            if not shared:
                continue
            algo = True
            tab = collections.Counter((cols[pa][r], cols[pb][r]) for r in shared)
            tot = sum(tab.values())
            print('posiciones %d y %d  |  %d lecturas abarcan ambas  |  '
                  'separacion %d pb' % (pa, pb, tot, abs(pb - pa)))
            for (x, y), n in tab.most_common():
                print('    %s-%s : %4d  (%.0f%%)' % (x, y, n, 100.0 * n / tot))
            dos = sum(n for _, n in tab.most_common(2))
            print('    las dos combinaciones mayoritarias suman %.0f%%' %
                  (100.0 * dos / tot))
            print('    -> %s' % ('LIGADAS, compatible con dos clones'
                                 if dos / tot >= 0.90 else
                                 'independientes, compatible con error de basecalling'))
            print()
    if not algo:
        print('Ninguna lectura abarca dos de esas posiciones: estan demasiado '
              'separadas para medir ligamiento con estas longitudes.')


if __name__ == '__main__':
    main()

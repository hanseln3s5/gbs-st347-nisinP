#!/usr/bin/env python3
"""
Reproduce las pruebas de asociacion del manuscrito.

Entrada: la tabla de un genoma por fila que produce generar_suplementarios.sh
(suplementarios/Tabla_S1_genomas.tsv), con las columnas:

    accession  ST  serotype  nisinP_operon  nsr  resistance_genes

Uso:
    python3 estadisticas.py suplementarios/Tabla_S1_genomas.tsv

Salida: cada tabla 2x2, su prueba exacta de Fisher a dos colas, y la razon de
momios con intervalo de confianza del 95%. Cuando una casilla es cero se aplica
la correccion de Haldane-Anscombe (se suma 0.5 a las cuatro casillas) antes de
calcular la razon de momios y su intervalo; la prueba de Fisher se hace sobre la
tabla sin corregir.
"""

import math
import sys
import collections

try:
    from scipy.stats import fisher_exact
except ImportError:
    sys.exit('Falta scipy. Activa el entorno gbs-asm, que lo tiene.')


def odds_ratio_ci(a, b, c, d, conf=0.95):
    """Razon de momios con IC, con correccion de Haldane-Anscombe si hace falta."""
    corrected = 0 in (a, b, c, d)
    if corrected:
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    orv = (a * d) / (b * c)
    se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    # z = 1.96 y no el valor exacto 1.959963984540054: es el que se uso al
    # calcular los intervalos publicados, y la diferencia (188,623 contra
    # 188,609 en el limite superior de tet(O)) solo aparece en intervalos tan
    # anchos como estos. Se deja asi para que el script reproduzca el
    # manuscrito al digito.
    z = 1.96
    lo = math.exp(math.log(orv) - z * se)
    hi = math.exp(math.log(orv) + z * se)
    return orv, lo, hi, corrected


def test(nombre, expuestos_si, expuestos_no, no_exp_si, no_exp_no):
    """expuestos = portadores del operon; si/no = presencia del gen."""
    a, b, c, d = expuestos_si, expuestos_no, no_exp_si, no_exp_no
    _, p = fisher_exact([[a, b], [c, d]])
    orv, lo, hi, corr = odds_ratio_ci(a, b, c, d)
    print('  %s' % nombre)
    print('                  presente  ausente')
    print('    portadores     %7d  %7d' % (a, b))
    print('    no portadores  %7d  %7d' % (c, d))
    print('    Fisher a dos colas  p = %.3g' % p)
    print('    OR = %.1f   IC95%% %.1f a %.1f%s'
          % (orv, lo, hi, '   (Haldane-Anscombe)' if corr else ''))
    print()


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    path = sys.argv[1]

    filas = []
    with open(path) as fh:
        head = fh.readline().rstrip('\n').split('\t')
        idx = {c: i for i, c in enumerate(head)}
        for col in ('accession', 'ST', 'nisinP_operon', 'nsr',
                    'resistance_genes'):
            if col not in idx:
                sys.exit('La tabla no tiene la columna %r. Columnas: %s'
                         % (col, head))
        for line in fh:
            if line.strip():
                filas.append(line.rstrip('\n').split('\t'))

    def g(r, c):
        return r[idx[c]]

    n = len(filas)
    print('genomas en la tabla: %d\n' % n)
    if n != 83:
        print('AVISO: el manuscrito analiza 83 genomas y aqui hay %d.\n' % n)

    portador = [r for r in filas if g(r, 'nisinP_operon') == 'yes']
    otros = [r for r in filas if g(r, 'nisinP_operon') != 'yes']
    print('portadores del operon de nisina P: %d   resto: %d\n'
          % (len(portador), len(otros)))

    def tiene(r, gen):
        return gen.lower() in g(r, 'resistance_genes').lower()

    print('=== asociacion entre el operon y los genes de resistencia ===\n')
    for gen, etiqueta in (('tet(O)', 'tet(O)'), ('erm(B)', 'erm(B)')):
        clave = gen.replace('(', '').replace(')', '')
        si_p = sum(1 for r in portador if tiene(r, clave))
        si_o = sum(1 for r in otros if tiene(r, clave))
        test(etiqueta, si_p, len(portador) - si_p, si_o, len(otros) - si_o)

    print('=== nsr por linaje ===\n')
    con_nsr = sum(1 for r in filas if g(r, 'nsr') == 'yes')
    print('  genomas con nsr: %d de %d' % (con_nsr, n))
    sin_nsr_st = collections.Counter(g(r, 'ST') for r in filas
                                     if g(r, 'nsr') != 'yes')
    print('  ST de los genomas SIN nsr: %s'
          % ', '.join('ST%s x%d' % (k, v) for k, v in
                      sorted(sin_nsr_st.items())))
    st1 = [r for r in filas if g(r, 'ST') == '1']
    if st1:
        st1_sin = sum(1 for r in st1 if g(r, 'nsr') != 'yes')
        st1_op = sum(1 for r in st1 if g(r, 'nisinP_operon') == 'yes')
        print('  ST1: %d genomas, %d sin nsr, %d con el operon'
              % (len(st1), st1_sin, st1_op))
    print()

    print('=== ST1010: determinante de aminoglucosidos ===\n')
    st1010 = [r for r in filas if g(r, 'ST') == '1010']
    if st1010:
        hlgr = sum(1 for r in st1010
                   if 'aac(6' in g(r, 'resistance_genes').lower()
                   or 'aph(2' in g(r, 'resistance_genes').lower())
        print('  ST1010: %d genomas, %d con aac(6\')-Ie/aph(2\'\')-Ia'
              % (len(st1010), hlgr))
    else:
        print('  no hay genomas ST1010 en la tabla')
    print()

    print('Las fracciones alelicas de las posiciones discordantes entre S1 y')
    print('S2 y su prueba de homogeneidad se calculan en fracciones.py.')


if __name__ == '__main__':
    main()

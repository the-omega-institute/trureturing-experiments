#!/usr/bin/env python3
"""Exact Definition 3.1 checks using Singular's characteristic-independent radical().

Python generates every tableau, support, and partition; Singular computes ideal
intersections, radicals, and two-way Groebner containment. No numerical evidence.
Requires Python 3 and Singular (primdec.lib). Run: python3 check.py
"""
from itertools import combinations, permutations
import re
import shutil
import subprocess
import tempfile
from pathlib import Path


def partitions(n, cap=None):
    if n == 0:
        yield ()
        return
    for k in range(min(n, n if cap is None else cap), 0, -1):
        for rest in partitions(n-k, k):
            yield (k,) + rest


def specht(shape, support, r):
    cells = [(i, j) for i, size in enumerate(shape) for j in range(size)]
    result = set()
    for labels in permutations(support):
        t = dict(zip(cells, labels))
        factors = [f'(x{t[i,j]}^{r}-x{t[k,l]}^{r})'
                   for i, j in cells for k, l in cells if j == l and i < k]
        result.add('*'.join(factors) or '1')
    return sorted(result)


def ideal_text(gens):
    return ','.join(gens) or '0'


def source_ideal(n, z, shape, r):
    m = n-z
    code = ['ideal I=1;']
    for support in combinations(range(1, n+1), m):
        gens = [f'x{i}' for i in range(1,n+1) if i not in support]
        gens += specht(shape, support, r)
        code += ['ideal J='+ideal_text(gens)+';', 'I=intersect(I,J);', 'kill J;']
    squarefree = ['*'.join(f'x{i}' for i in s) or '1'
                  for s in combinations(range(1,n+1), m)]
    code += ['ideal Q='+ideal_text(squarefree)+';', 'I=std(intersect(I,Q));']
    return '\n'.join(code)


def run_case(n, z, shape, r, char):
    coeff = '(0,a)' if char == 0 and r == 3 else str(char)
    variables = ','.join(f'x{i}' for i in range(1,n+1)) or 'dummy'
    code = f'LIB "primdec.lib";\nring R={coeff},({variables}),dp;\n'
    if coeff == '(0,a)':
        code += 'minpoly=a^2+a+1;\n'
    code += source_ideal(n,z,shape,r)
    code += '\nideal RI=radical(I);\nprint("RAD="+string(idealsEqual(I,RI)));\n'
    if n == 3 and z == 1 and shape == (1,1):
        expected = ','.join(f'x{i}*x{j}*(x{i}^{r}-x{j}^{r})'
                            for i,j in combinations(range(1,4),2))+',x1*x2*x3'
        code += f'ideal E={expected};\nprint("FID="+string(idealsEqual(I,E)));\n'
        if char == r:
            # Global witness, not merely a localized or pointwise test.
            code += ('poly w=x1*x2*(x1-x2);\n'
                     'print("W_NONMEM="+string(reduce(w,I)!=0));\n'
                     f'print("W_POWER_MEM="+string(reduce(w^{r},I)==0));\n')
    code += 'quit;\n'
    with tempfile.TemporaryDirectory(prefix='singular-', dir=Path(__file__).parent) as td:
        program = Path(td)/'case.sing'; program.write_text(code)
        p = subprocess.run([shutil.which('Singular'), '-q', str(program)],
                           capture_output=True, text=True, timeout=180)
    assert p.returncode == 0, (p.returncode,p.stdout,p.stderr)
    assert '?' not in p.stdout and 'error' not in p.stdout.lower(), p.stdout
    readings = dict(re.findall(r'(RAD|FID|W_NONMEM|W_POWER_MEM)=(\d+)',p.stdout))
    expected_radical = char == 0 or len(shape) <= 1
    assert readings.get('RAD') == str(int(expected_radical)), (n,z,shape,r,char,p.stdout)
    if 'FID' in readings:
        assert readings['FID'] == '1', p.stdout
    if char == r and n == 3 and z == 1 and shape == (1,1):
        assert readings.get('W_NONMEM') == readings.get('W_POWER_MEM') == '1',p.stdout
    field = ('Q(zeta_3)' if r == 3 else 'Q=Q(zeta_2)') if char == 0 else f'F_{char}'
    suffix = ' fidelity=PASS' if 'FID' in readings else ''
    if 'W_NONMEM' in readings:
        suffix += f' witness=x1*x2*(x1-x2) nonmember=1 power_{r}_member=1'
    print(f'field={field} r={r} n={n} z={z} mu={shape} radical={readings["RAD"]}{suffix}',flush=True)
    return expected_radical


def main():
    assert shutil.which('Singular'), 'Singular is required'
    print(subprocess.check_output(['Singular','--version'],text=True).splitlines()[0])
    print('METHOD: primdec.lib radical(I), equality by two-way exact Groebner containment')
    print('RANGE: r=2,3; n=0,1,2,3,4; z=0,...,n; every partition of n-z')
    total = radical_count = 0
    for r in (2,3):
        for char in (0,r):
            for n in range(5):
                for z in range(n+1):
                    for shape in partitions(n-z):
                        radical_count += run_case(n,z,shape,r,char)
                        total += 1
    print(f'TOTAL cases={total} radical={radical_count} nonradical={total-radical_count}')
    print('FIDELITY: source Examples 4.1 and 4.3 polynomial ideal presentations reproduced')
    print('BOUNDARY: finite fields test the literal polynomial definition, not the standing infinite-field group realization')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Exact finite checks for Faber--Schober, arXiv:2609.01501v1, Question 4.5.

Requires Python >= 3.12 and Singular (primdec.lib, minAssGTZ). No floating point.
The geometric ideal is formed from the source's equations, before decomposition.
Coordinate permutation identifies each ordered tuple with its sorted tuple;
all ordered tuples are checked, and every sorted orbit is decomposed over Q/F2.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, combinations_with_replacement, product
from math import comb, prod
from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile


def components(q, j):
    out = []
    for active in combinations(range(len(q)), j):
        a = frozenset(active)
        unit = [i for i in active if q[i] == 1]
        for axes in product((0, 1), repeat=len(unit)):
            out.append((a, tuple(zip(unit, axes))))
    return out


def codim(q, component):
    active, _ = component
    return sum(2*d for i, d in enumerate(q) if i not in active) + len(active)


def model(q, j):
    return Counter(codim(q, c) for c in components(q, j))


def enumerator(q):
    p = {(0, 0): 1}
    for d in q:
        nxt = defaultdict(int)
        for (j, c), v in p.items():
            nxt[j, c + 2*d] += v
            nxt[j + 1, c + 1] += v * (2 if d == 1 else 1)
        p = nxt
    return p


def multiplicity_formula(q, j):
    m = Counter(q)
    ds = sorted(m)
    out = Counter()
    for ks in product(*(range(m[d] + 1) for d in ds)):
        if sum(ks) != len(q) - j:
            continue
        k = dict(zip(ds, ks))
        c = 2*sum(d*k[d] for d in ds) + j
        out[c] += 2**(m[1] - k.get(1, 0)) * prod(comb(m[d], k[d]) for d in ds)
    return out


def source_count(q, j):
    n, alpha = len(q), q.count(1)
    def choose(a, b):
        return comb(a, b) if 0 <= b <= a else 0
    return sum(choose(n-alpha, n-j-k)*choose(alpha, k)*2**(alpha-k)
               for k in range(alpha + 1))


def source_total(q):
    n, alpha = len(q), q.count(1)
    return 3**alpha * 2**(n-alpha) - (n+1)*2**alpha + alpha*Fraction(2)**(alpha-1)


def contained(c, d):
    # Variety inclusion: active blocks only increase; old unit axes must agree.
    a, ax = c
    b, bx = d
    return a <= b and all(dict(bx).get(i) == v for i, v in ax)


def union_pairs(q):
    n = len(q)
    out = Counter()
    for zero in combinations(range(n), 2):
        c = n-2 + 2*sum(q[i] for i in zero)
        out[c] += prod(2 if q[i] == 1 else 1 for i in range(n) if i not in zero)
    return out


def singular_check(q, char, executable, scratch):
    n = len(q)
    blocks = [[f'z{a}_{i}' for i in range(2*d)] for a, d in enumerate(q)]
    names = [v for block in blocks for v in block]
    gs = ['+'.join(f'{block[2*i]}*{block[2*i+1]}' for i in range(q[a]))
          for a, block in enumerate(blocks)]
    diffs = [f'({gs[0]})-({gs[a]})' for a in range(1, n)]
    lines = ['LIB "primdec.lib";', f'ring r={char},({",".join(names)}),dp;',
             'ideal I,K,J,U,E;', 'list P;', 'int j,h,k,found,ok;', 'ok=1;']
    for j in range(n-1):
        zero_sets = list(combinations(range(n), n-j))
        for h, zero in enumerate(zero_sets):
            generators = [v for a in zero for v in blocks[a]] + diffs
            lines.append('K='+','.join(generators)+';')
            lines.append('I=K;' if h == 0 else 'I=intersect(I,K);')
        lines.append('I=std(I);')
        if j == 0:
            lines.append('U=I;')
        else:
            # The source's literal union is assembled independently by intersection.
            lines.append('U=std(intersect(U,I));')
        lines += ['P=minAssGTZ(I);',
                  f'if(size(P)!={len(components(q,j))}){{ok=0;}}']
        for component in components(q, j):
            active, axes = component
            gens = [v for a in range(n) if a not in active for v in blocks[a]]
            for a in sorted(active):
                if q[a] == 1:
                    gens.append(blocks[a][dict(axes)[a]])
                else:
                    gens.append(gs[a])
            lines += ['E=std(ideal('+','.join(gens)+'));', 'found=0;',
                      'for(k=1;k<=size(P);k++){', 'J=std(P[k]);',
                      'if(size(reduce(E,J))==0 && size(reduce(J,E))==0){',
                      'found=found+1;',
                      f'if(nvars(r)-dim(J)!={codim(q, component)}){{ok=0;}}',
                      '}', '}', 'if(found!=1){ok=0;}']
        lines.append(f'if(ok==0){{print("CHECK_FAILED stage {j}");}}')
    lines += ['if(size(reduce(U,I))!=0 || size(reduce(I,U))!=0){ok=0;}',
              'if(ok==1){print("SINGULAR_OK");}', 'quit;']
    script = scratch / 'case.sing'
    script.write_text('\n'.join(lines)+'\n')
    run = subprocess.run([executable, '-q', str(script)], text=True, capture_output=True, timeout=240)
    if run.returncode or run.stderr.strip() or run.stdout.strip() != 'SINGULAR_OK':
        raise AssertionError((q, char, run.returncode, run.stdout, run.stderr))
    return n-1


def main():
    assert sys.version_info >= (3, 12)
    executable = shutil.which('Singular')
    assert executable, 'Singular required: finite counting alone is not a geometry check'
    print('Python', sys.version.split()[0])
    version = subprocess.run([executable, '--version'], input='quit;\n', text=True,
                             capture_output=True, check=True).stdout.splitlines()[0]
    print(version)
    print('Source: Faber--Schober, arXiv:2609.01501v1, Proposition 4.3, Corollary 4.4, Question 4.5')
    ordered = stages = 0
    for n in range(2, 5):
        for q in product(range(1, 4), repeat=n):
            ordered += 1
            p = enumerator(q)
            total = 0
            # Verify the actual coordinate permutation used for sorted-orbit reduction.
            permutation = sorted(range(n), key=lambda a: (q[a], a))
            sorted_q = tuple(q[a] for a in permutation)
            to_sorted = {a: b for b, a in enumerate(permutation)}
            maximal = components(q, n-2)
            for j in range(n-1):
                stages += 1
                cs = components(q, j)
                permuted = {(frozenset(to_sorted[a] for a in active),
                             tuple(sorted((to_sorted[a], side) for a, side in axes)))
                            for active, axes in cs}
                assert permuted == set(components(sorted_q, j))
                expected = model(q, j)
                coefficient = Counter({c:v for (jj,c),v in p.items() if jj == j})
                assert expected == coefficient == multiplicity_formula(q, j)
                assert sum(expected.values()) == source_count(q, j)
                assert all(any(contained(c, d) for d in maximal) for c in cs)
                assert not any(c != d and contained(c, d) for c in cs for d in cs)
                if len(set(q)) == 1:
                    d = q[0]
                    assert set(expected) == {2*d*(n-j) + j}
                    assert sum(expected.values()) == comb(n,j)*(2**j if d == 1 else 1)
                total += sum(expected.values())
            assert total == source_total(q)
            assert model(q, n-2) == union_pairs(q)
            # Stage-tagged polynomial: delete terms with zero/one zero blocks.
            stage_poly = Counter()
            for (j,c),v in p.items():
                if j <= n-2:
                    stage_poly[c] += v
            direct = sum((model(q,j) for j in range(n-1)), Counter())
            assert stage_poly == direct
    print(f'Ordered tuples: {ordered}; stage distributions: {stages}; exact model/coefficient/support/fidelity checks passed')
    algebra_stages = orbits = 0
    with tempfile.TemporaryDirectory(prefix='fs2609q45-') as folder:
        scratch = Path(folder)
        for n in range(2, 5):
            for q in combinations_with_replacement(range(1,4), n):
                orbits += 1
                for char in (0, 2):
                    algebra_stages += singular_check(q, char, executable, scratch)
                print('Minimal-prime and ambient-codimension checks Q/F2:',q, flush=True)
    print(f'Permutation orbits: {orbits}; Singular stage decompositions: {algebra_stages}')
    example = union_pairs((1,1,4,4))
    assert example == Counter({6:1,12:8,18:4})
    print('Support-gap illustration q=(1,1,4,4):',json.dumps(dict(sorted(example.items()))))
    print('Scope: finite exact checks, not a universal proof or Lean verification')
    print('check.py sha256:',hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

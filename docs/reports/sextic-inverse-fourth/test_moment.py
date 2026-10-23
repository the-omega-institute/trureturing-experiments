#!/usr/bin/env python3
"""Independent scalar/full-unit expansion against the batched moment engine."""
import cmath
from collections import Counter, defaultdict
import math
import sys
sys.dont_write_bytecode = True
from eisenstein import Engine, UNITS, mul, norm
from moment import run_scale


def weight(x, c, r):
    t = (x - c) / r
    return math.exp(1 - 1 / (1 - t*t)) if abs(t) < 1 else 0.0


def value(engine, d, u):
    e = engine.symbol_exponent(d, u)
    return 0j if e < 0 else cmath.exp(math.pi * 1j * e / 3)


def reference(engine, U):
    D = math.sqrt(U)
    columns = [d for d in engine.ideals if D/2 < norm(d) < 1.5*D]
    weights = {d: weight(norm(d)/D, 1, .5) for d in columns}
    g = defaultdict(float)
    core = {}
    head = set()
    active = [d for d in columns if engine.mu(d)]
    for d in active:
        for e in active:
            a = mul(d, e)
            g[a] += engine.mu(d)*engine.mu(e)*weights[d]*weights[e]/D
            f = Counter(engine.factors[d]) + Counter(engine.factors[e])
            core[a] = tuple(sorted(p for p, v in f.items() if v % 2))
            if not set(engine.factors[d]).intersection(engine.factors[e]):
                head.add(a)
    total = defaultdict(float)
    for p in engine.ideals:
        if not U < norm(p) < 2*U or not engine.is_sixth_free(p):
            continue
        w = weight(norm(p)/U, 1.5, .5)
        for unit in UNITS:
            u = mul(p, unit)
            total['row_mass'] += w
            character = {d: value(engine, d, u) for d in columns}
            for name, coeff in [
                ('mobius', engine.mu),
                ('liouville', engine.liouville),
                ('ordinary', lambda d: 1),
                ('squarefree_positive', lambda d: abs(engine.mu(d))),
            ]:
                z = sum(coeff(d)*weights[d]*character[d] for d in columns)/math.sqrt(D)
                total[name+'_second'] += w*abs(z)**2
                total[name+'_fourth'] += w*abs(z)**4
            by_core = defaultdict(complex)
            zhead = 0j
            zg = 0j
            diag = 0.0
            headdiag = 0.0
            for a, c in g.items():
                z = c*value(engine, a, u)
                by_core[core[a]] += z
                zg += z
                diag += abs(z)**2
                if a in head:
                    zhead += z
                    headdiag += abs(z)**2
            total['h2'] += w*sum(abs(z)**2 for z in by_core.values())
            total['diagonal'] += w*diag
            total['t1_diagonal'] += w*headdiag
            total['t1_energy'] += w*abs(zhead)**2
            total['t_gt1_energy'] += w*abs(zg-zhead)**2
    return total


def main():
    engine = Engine(2048)
    cases = 0
    for U in [256, 1024]:
        expected = reference(engine, U)
        for chunk in [37, 256]:
            got = run_scale(engine, U, chunk)
            for key, wanted in expected.items():
                error = abs(got[key] - wanted)
                if error > 2e-10 * max(1, abs(wanted)):
                    raise AssertionError((U, chunk, key, got[key], wanted))
                cases += 1
    print(f'{cases} scalar/full-unit moment comparisons passed')


if __name__ == '__main__':
    main()

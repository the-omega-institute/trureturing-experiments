#!/usr/bin/env python3
"""Exact, finite checks for Hanna's OEIS A266481 prime-multiple conjecture.

Sources: https://oeis.org/A266481 (revision 35),
https://oeis.org/A266526 (revision 9), https://oeis.org/A266521 (revision 13).
The checked finite statements do not establish the infinite conjecture.
Python 3 standard library only; every arithmetic operation is exact.
"""
from fractions import Fraction as Q
from itertools import combinations, product
from math import comb, factorial

BOUND = 30
SOURCE_A = [
    1, 1, 5, 55, 993, 25501, 857773, 35850795, 1795564865,
    104972371417, 7022842421301, 529428563641759, 44421725002096225,
    4106744812439019765, 414834196219620026333, 45462732300569936279251,
    5373006006732947705188737, 681229881246574750274962225,
    92237589983019368975021777125, 13283769418970268811752725081607,
    2027649185923009220298941142143201,
    326999803592314489529958494308640461,
    55558592280735155060861740192416874125,
]
SOURCE_T = [
    1, 4, 42, 752, 19360, 654912, 27546736, 1388207872, 81621893376,
    5488951731200, 415721105434624, 35026876903256064,
    3250356630453317632, 329437813126362185728, 36214170617862339840000,
    4291812357982293898231808, 545518054282041342531076096,
    74032137722410904128877494272, 10684317262536125210489796296704,
    1634019721630446295055397683200000,
]
SOURCE_C = [
    [1, 2, 1], [15, 28, 18, 4], [683, 1278, 933, 316, 42],
    [62038, 117440, 92680, 38240, 8272, 752],
    [9342629, 17880090, 14855385, 6881640, 1880340, 288048, 19360],
    [2100483216, 4054752672, 3490688496, 1743156480, 547098240,
     108228192, 12523584, 654912],
]


def mul(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i+j] += x*y
    return c


def source_connected_polynomials(bound):
    """C_m=m![x^m]log(sum (m+y)^(2m)x^m/m!), by F'=C'F."""
    f = [[comb(2*m, j) * m**(2*m-j) for j in range(2*m+1)]
         for m in range(bound+1)]
    c = [[]]
    for m in range(1, bound+1):
        row = f[m][:]
        for k in range(1, m):
            factor = comb(m-1, k-1)
            for j, z in enumerate(mul(c[k], f[m-k])):
                row[j] -= factor*z
        while row and row[-1] == 0:
            row.pop()
        assert len(row) <= m+2, (m, row)
        assert all(x >= 0 for x in row), m
        c.append(row)
    return c


def assemble(t):
    """Exponential formula recurrence for labelled SET of connected components."""
    a = [1]
    for m in range(1, len(t)):
        a.append(sum(comb(m-1, k-1)*t[k]*a[m-k] for k in range(1, m+1)))
    return a


def ladd(a, b):
    c = a.copy()
    for e, v in b.items():
        c[e] = c.get(e, Q(0)) + v
        if not c[e]:
            del c[e]
    return c


def lmul(a, b):
    c = {}
    for e, x in a.items():
        for f, y in b.items():
            c[e+f] = c.get(e+f, Q(0)) + x*y
    return {e: v for e, v in c.items() if v}


def direct_source_limit(bound):
    """Direct NAME expansion: t=1/N, (1+U)^t=sum choose(t,k)U^k.

    Its coefficients are Laurent polynomials in t. Negative powers must cancel;
    the coefficient of t^0 is the exact coefficientwise limit as t tends to 0.
    This does not use the connected-polynomial degree bound.
    """
    u = [{}] + [
        {j-m: Q(comb(2*m, j)*m**j, factorial(m)) for j in range(2*m+1)}
        for m in range(1, bound+1)]
    power = [{0: Q(1)}] + [{} for _ in range(bound)]
    choose = {0: Q(1)}
    result = [{} for _ in range(bound+1)]
    for k in range(bound+1):
        for m in range(bound+1):
            result[m] = ladd(result[m], lmul(choose, power[m]))
        nxt = [{} for _ in range(bound+1)]
        for m in range(1, bound+1):
            for j in range(1, m+1):
                nxt[m] = ladd(nxt[m], lmul(power[m-j], u[j]))
        power = nxt
        choose = {e: v/Q(k+1) for e, v in
                  lmul(choose, {1: Q(1), 0: Q(-k)}).items() if v}
    answer = []
    for m, row in enumerate(result):
        assert all(e >= 0 for e in row), (m, row)
        value = row.get(0, Q(0))*factorial(m)
        assert value.denominator == 1, (m, value)
        answer.append(value.numerator)
    return answer


def prufer_tree(n, code):
    degree = [1] * n
    for v in code:
        degree[v] += 1
    edges = []
    for v in code:
        leaf = next(i for i, d in enumerate(degree) if d == 1)
        edges.append((leaf, v))
        degree[leaf] -= 1
        degree[v] -= 1
    leaves = [i for i, d in enumerate(degree) if d == 1]
    edges.append(tuple(leaves))
    return edges


def orientation_weight(n, edges, directions):
    outdegree = [0] * n
    for (u, v), direction in zip(edges, directions):
        outdegree[v if direction else u] += 1
    # Injections from outgoing edges to two distinguished slots: (2)_d.
    factors = [1, 2, 2]
    if max(outdegree, default=0) > 2:
        return 0
    weight = 1
    for d in outdegree:
        weight *= factors[d]
    return weight


def tree_count(n):
    if n == 1:
        return 1
    total = 0
    for code in product(range(n), repeat=n-2):
        edges = prufer_tree(n, code)
        total += sum(orientation_weight(n, edges, ds)
                     for ds in product(range(2), repeat=n-1))
    return total


def undirected_forests(n):
    pairs = list(combinations(range(n), 2))
    for mask in range(1 << len(pairs)):
        if mask.bit_count() >= n and n:
            continue
        edges = [e for i, e in enumerate(pairs) if (mask >> i) & 1]
        parents = list(range(n))
        def root(v):
            while parents[v] != v:
                v = parents[v]
            return v
        for u, v in edges:
            ru, rv = root(u), root(v)
            if ru == rv:
                break
            parents[ru] = rv
        else:
            yield edges


def slot_forests(n):
    """Enumerate actual maps V x {0,1}->Option(V), not just weighted graphs."""
    for edges in undirected_forests(n):
        for ds in product(range(2), repeat=len(edges)):
            outgoing = [[] for _ in range(n)]
            for (u, v), d in zip(edges, ds):
                source, target = (v, u) if d else (u, v)
                outgoing[source].append(target)
            if any(len(ts) > 2 for ts in outgoing):
                continue
            options = []
            for ts in outgoing:
                if not ts:
                    options.append([(-1, -1)])
                elif len(ts) == 1:
                    options.append([(ts[0], -1), (-1, ts[0])])
                else:
                    options.append([(ts[0], ts[1]), (ts[1], ts[0])])
            for choices in product(*options):
                yield tuple(t for pair in choices for t in pair)


def relabel(f, permutation):
    image = [-1] * len(f)
    for v in range(len(permutation)):
        for s in range(2):
            target = f[2*v+s]
            image[2*permutation[v]+s] = -1 if target == -1 else permutation[target]
    return tuple(image)


def main():
    c = source_connected_polynomials(BOUND)
    t = [0] + [row[m+1] if len(row) > m+1 else 0
               for m, row in enumerate(c[1:], 1)]
    a = assemble(t)
    assert c[1:7] == SOURCE_C
    assert t[1:len(SOURCE_T)+1] == SOURCE_T
    assert a[:len(SOURCE_A)] == SOURCE_A
    print('FIDELITY A266481: a(0)..a(22) match all 23 displayed terms')
    print('FIDELITY A266526: T(1)..T(20) match all 20 displayed terms')
    print('FIDELITY A266521: connected rows 1..6 match exactly')
    print('SOURCE_NOTE A266521 example exp(A(x,y)) prints (1+y)*x; '
          'its NAME and first connected row instead give (1+y)^2*x')
    direct = direct_source_limit(6)
    assert direct == a[:7]
    print('DIRECT_NAME_LIMIT a(0)..a(6):', direct)
    print('CONNECTED_DEGREE_BOUND checked for 1<=m<=30')
    trees = [tree_count(m) for m in range(1, 7)]
    assert trees == t[1:7]
    print('PRUFER_SLOT_TREES T(1)..T(6):', trees)
    forests = {}
    for m in range(6):
        fs = list(slot_forests(m))
        assert len(fs) == len(set(fs)) == a[m], m
        weighted = sum(orientation_weight(m, edges, ds)
                       for edges in undirected_forests(m)
                       for ds in product(range(2), repeat=len(edges)))
        assert weighted == len(fs), m
        forests[m] = fs
    print('ALL_SLOT_FORESTS a(0)..a(5):', [len(forests[m]) for m in range(6)])
    for p, blocks in [(2, 1), (2, 2), (3, 1), (5, 1)]:
        m = p*blocks
        generators = []
        for i in range(blocks):
            perm = list(range(m))
            for j in range(p):
                perm[p*i+j] = p*i+(j+1) % p
            generators.append(perm)
        fixed = [f for f in forests[m]
                 if all(relabel(f, g) == f for g in generators)]
        assert fixed == [(-1,)*(2*m)], (p, blocks, fixed)
        # Also check closure under each group generator, independently of counting.
        universe = set(forests[m])
        assert all(relabel(f, g) in universe for f in universe for g in generators)
        print(f'GLOBAL_FIXED_FORESTS p={p}, blocks={blocks}, total={len(universe)}, fixed=1')
    assert (1, -1, 0, -1) not in set(forests[2])
    print('NEGATIVE_CONTROL reciprocal two-vertex arrows are rejected')
    primes = [p for p in range(2, BOUND+1)
              if all(p % d for d in range(2, p))]
    checks = [(p, n) for p in primes for n in range(BOUND//p+1)]
    assert all((a[p*n]-1) % p == 0 for p, n in checks)
    print(f'PRIME_CONGRUENCE {len(checks)} pairs: every prime p<=30 and n>=0 with p*n<=30')
    assert Q(t[4], factorial(3)) == Q(376, 3)
    print('PRIOR_ART_HYPOTHESIS logarithmic weight w_4 = 376/3 is not integral')
    print('Finite checks only; the all-index theorem remains a preregistered proof obligation.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

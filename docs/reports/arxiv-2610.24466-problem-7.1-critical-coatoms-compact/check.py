#!/usr/bin/env python3
"""Exact finite checks for Wang–Dou arXiv:2610.24466v1 Problem 7.1.

This independently implements the available-profile semantics and Table 1.
It checks the proposed finite centrality/guard basis on finite factors and
products, and exhausts all MCL expansions of the four-element Boolean algebra.
It does not decide compactness of an infinite logic lattice.
"""
from itertools import product

IA_TABLE = {(0, 1): (1,), (1, 0): (0,), (1, 1): (1,),
            (1, 2): (0,), (2, 1): (0,), (2, 2): (1,)}
TABLES = {
    'IA': IA_TABLE,
    'Det': {(0, 0): (1,), (1, 0): (1,), (0, 1): (1,), (1, 1): (0, 2)},
    'Mix': dict(IA_TABLE, **{}) | {(3, 1): (1,), (3, 3): (2, 3)},
}
ACTION_COUNTS = {'IA': 3, 'Det': 2, 'Mix': 4}
MIN_PRIME = {'IA': 2, 'Det': 3, 'Mix': 5}


def subset(x, y):
    return x & y == x


def shift(x, p, k=1):
    return sum(1 << s for s in range(p) if x & (1 << ((s + k) % p)))


def table_op(kind, c, x, p):
    g1, g2, g3 = (shift(x, p, j) for j in (1, 2, 3))
    active = c & 3
    if active == 0:
        return {'IA': x & g1, 'Det': x & g1 & g2,
                'Mix': x & g1 & g2 & g3}[kind]
    if active == 1:
        return g1
    if active == 2:
        return {'IA': x, 'Det': g1, 'Mix': x | (g2 & g3)}[kind]
    return {'IA': x | g1, 'Det': g1 | (x & g2),
            'Mix': x | g1 | (g2 & g3)}[kind]


def semantic_operations(kind, n, p):
    """Full common action alphabet, with every nonzero dummy action unavailable."""
    size = 1 << p
    operations = [[0] * size for _ in range(1 << n)]
    for s in range(p):
        available = []
        for pi in product(range(ACTION_COUNTS[kind]), repeat=n):
            offsets = TABLES[kind].get(pi[:2], ()) if all(v == 0 for v in pi[2:]) else ()
            outcome = sum(1 << t for t in {(s + j) % p for j in offsets})
            if outcome:
                available.append((pi, outcome))
        for c in range(1 << n):
            aggregates = {}
            for pi, outcome in available:
                restriction = tuple(pi[i] for i in range(n) if c & (1 << i))
                aggregates[restriction] = aggregates.get(restriction, 0) | outcome
            for x in range(size):
                if any(subset(outcome, x) for outcome in aggregates.values()):
                    operations[c][x] |= 1 << s
    return operations


def q(x, operations, p, top):
    result = top
    for _ in range(p):
        result &= x
        x = operations[1][x]
    return result


def agreement(x, y, top):
    """Boolean truth value of x approximately equals y."""
    return top ^ (x ^ y)


def delta_values(operations, top):
    """All truth values for each Ser, IA, Det equation, separately."""
    nc, size = len(operations), top + 1
    result = []
    for c in range(nc):
        result.append({operations[c][top]})
    for c in range(nc):
        for d in range(nc):
            if c & d == 0:
                result.append({(top ^ (operations[c][x] & operations[d][y])) |
                               operations[c | d][x & y]
                               for x in range(size) for y in range(size)})
    for c in range(nc):
        result.append({(top ^ operations[c][x | y]) | operations[c][x] | operations[nc-1][y]
                       for x in range(size) for y in range(size)})
    return result


def epsilon_values(operations, top, p, kind):
    """Finite basis Delta_A: Boolean shift automorphism, period, Table 1."""
    size, g = top + 1, operations[1]
    result = [{agreement(g[x & y], g[x] & g[y], top)
               for x in range(size) for y in range(size)},
              {agreement(g[top ^ x], top ^ g[x], top) for x in range(size)}]
    period = []
    for x in range(size):
        y = x
        for _ in range(p):
            y = g[y]
        period.append(agreement(y, x, top))
    result.append(set(period))
    for c in range(len(operations)):
        values = set()
        for x in range(size):
            g1, g2, g3 = g[x], g[g[x]], g[g[g[x]]]
            active = c & 3
            if active == 0:
                rhs = {'IA': x & g1, 'Det': x & g1 & g2, 'Mix': x & g1 & g2 & g3}[kind]
            elif active == 1:
                rhs = g1
            elif active == 2:
                rhs = {'IA': x, 'Det': g1, 'Mix': x | (g2 & g3)}[kind]
            else:
                rhs = {'IA': x | g1, 'Det': g1 | (x & g2),
                       'Mix': x | g1 | (g2 & g3)}[kind]
            values.add(agreement(operations[c][x], rhs, top))
        result.append(values)
    return result


def check_E(operations, top, p, deltas, epsilons):
    """All valuations, reduced to finite image sets; guard tuples are disjoint."""
    c_values = {q(v, operations, p, top) for values in deltas for v in values}
    for c in c_values:
        for eta in (c, top ^ c):
            for f in operations:
                for z in range(top + 1):
                    if f[z] & eta != f[z & eta] & eta:
                        return False
        for values in epsilons:
            for eps in values:
                if c | eps != top:
                    return False
    return True


def product_operations(left, right, left_bits):
    return [[f[x & ((1 << left_bits)-1)] | (g[x >> left_bits] << left_bits)
             for x in range(len(f)*len(g))] for f, g in zip(left, right)]


def product_values(left, right, left_bits):
    return [{x | (y << left_bits) for x in a for y in b} for a, b in zip(left, right)]


def is_mcl(operations, top):
    nc, size = len(operations), top + 1
    for c, f in enumerate(operations):
        if f[0] != 0 or any(not subset(v, operations[0][top]) for v in f):
            return False
        for d in range(nc):
            if subset(c, d) and any(not subset(f[x], operations[d][x]) for x in range(size)):
                return False
        for x in range(size):
            for y in range(size):
                if not subset(operations[0][(top ^ x) | y] & f[x], f[y]):
                    return False
    return True


def finite_checks():
    models = table_entries = scheme_valuations = products_checked = 0
    for kind in ('IA', 'Det', 'Mix'):
        for p in (2, 3, 5, 7):
            if p < MIN_PRIME[kind]:
                continue
            top = (1 << p)-1
            for n in (2, 3, 4):
                operations = semantic_operations(kind, n, p)
                nc = 1 << n
                for c in range(nc):
                    for x in range(top + 1):
                        assert operations[c][x] == table_op(kind, c, x, p), (kind, n, p, c, x)
                        table_entries += 1
                assert all(q(x, operations, p, top) == (top if x == top else 0)
                           for x in range(top + 1))
                deltas = delta_values(operations, top)
                epsilons = epsilon_values(operations, top, p, kind)
                assert len(deltas) == 2**(n+1) + 3**n
                assert len(epsilons) == 2**n + 3
                assert all(values == {top} for values in epsilons)
                ser, ia, det = deltas[:nc], deltas[nc:-nc], deltas[-nc:]
                valid = lambda equations: all(values == {top} for values in equations)
                assert valid(ser)
                assert valid(ia) == (kind == 'Det')
                assert valid(det) == (kind == 'IA')
                assert check_E(operations, top, p, deltas, epsilons)
                # A one-atom classical factor: every coalition operation is identity.
                classical = [[0, 1] for _ in range(nc)]
                cd = delta_values(classical, 1)
                ce = epsilon_values(classical, 1, p, kind)
                mixed = product_operations(operations, classical, p)
                md = product_values(deltas, cd, p)
                me = product_values(epsilons, ce, p)
                assert check_E(mixed, 2**(p+1)-1, p, md, me)
                # Two source factors, retaining the full coalition signature.
                source_squared = product_operations(operations, operations, p)
                sd = product_values(deltas, deltas, p)
                se = product_values(epsilons, epsilons, p)
                assert check_E(source_squared, 2**(2*p)-1, p, sd, se)
                products_checked += 2
                models += 1
                scheme_valuations += nc + (3**n + nc) * (top+1)**2
                print(f'finite family={kind} agents={n} prime={p}: Table1, type, Delta_A, E, two products OK')
    print(f'FIDELITY source_models={models} table_entries={table_entries} classical_axiom_valuations={scheme_valuations}')
    print(f'PROPOSED_BASIS factor_models={models} direct_products={products_checked}: all centrality and guards OK')


def exhaustive_four_element_check():
    """Exhaust all possible MCL operations on P({0,1}), for N={a,b}.

    MCL entails f_C(0)=0, f_C(1)=e=f_empty(1), f_C(x)<=e.
    Enumerate every tuple with these necessary properties, then filter by
    the four defining MCL equation families. For p=2, type IA, test E
    against CL or Delta_A on every resulting algebra. These are exactly
    the possible join members on this carrier: any nontrivial proper
    Boolean factor has two elements, and both source varieties have only
    classical two-element members. This finite equality is not a proof
    of the uniform result or of the general subdirect theorem.
    """
    candidates = mcl_models = e_models = classical_models = delta_models = 0
    for e in range(4):
        allowable = [x for x in range(4) if subset(x, e)]
        for entries in product(allowable, repeat=8):
            candidates += 1
            operations = [[0, entries[2*c], entries[2*c+1], e] for c in range(4)]
            if not is_mcl(operations, 3):
                continue
            mcl_models += 1
            deltas = delta_values(operations, 3)
            epsilons = epsilon_values(operations, 3, 2, 'IA')
            classical = all(v == {3} for v in deltas)
            delta = all(v == {3} for v in epsilons)
            holds = check_E(operations, 3, 2, deltas, epsilons)
            assert holds == (classical or delta), (operations, holds, classical, delta)
            e_models += int(holds)
            classical_models += int(classical)
            delta_models += int(delta)
    print(f'EXHAUSTIVE Boolean_carrier=4 agents=2 prime=2 type=IA candidates={candidates} MCL_models={mcl_models} E_models={e_models} CL_models={classical_models} Delta_models={delta_models}: E iff CL or Delta_A OK')


if __name__ == '__main__':
    finite_checks()
    exhaustive_four_element_check()
    print('ALL_CHECKS_OK')

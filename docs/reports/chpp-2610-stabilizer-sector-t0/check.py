#!/usr/bin/env python3
"""Exact finite checks for arXiv:2610.33221v1, Conjecture 5.6.

All configurations below extend to the infinite lattice as explicitly specified.
Only generators affected by a finite flip are counted. No cutoff total energy,
floating-point arithmetic, or inference of an infinite theorem from enumeration
is used. The infinite argument belongs to the preregistration prose.
"""
from itertools import product, combinations


def subsets(xs):
    xs = tuple(xs)
    for mask in range(1 << len(xs)):
        yield frozenset(x for i, x in enumerate(xs) if (mask >> i) & 1)


def parity(xs):
    return sum(xs) % 2


def energy_change(spin, flip, incident):
    """Eq. (4.7), by two finite violation counts (bits encode spin signs)."""
    nearby = set().union(*(incident(v) for v in flip)) if flip else set()
    affected = {g for g in nearby if len(g & flip) % 2}
    before = {g for g in affected if parity(spin(v) for v in g)}
    after = {g for g in affected if parity(spin(v) ^ (v in flip) for v in g)}
    return len(after) - len(before), affected


def line_incident(i):
    out = {frozenset((i - 1, i)), frozenset((i, i + 1))}
    if i == 0:
        out.add(frozenset((0,)))
    return out


def fidelity_biased_ising():
    """Prop. 6.23: two uniform configurations and one wall with satisfied pin.

    Enumerate every configuration on [-3,3] with each possible pair of constant
    infinite tails; test every finite spin flip supported in [-4,4]. The source
    criterion (at most one violated generator) is checked against this exhaustive
    finite test, and the four asymptotic/orbit labels are reproduced.
    """
    flips = tuple(subsets(range(-4, 5)))
    tested = boundary = 0
    sectors = set()
    for left, right in product((0, 1), repeat=2):
        for bits in product((0, 1), repeat=7):
            def spin(i):
                return left if i < -3 else right if i > 3 else bits[i + 3]
            walls = sum(spin(i) != spin(i + 1) for i in range(-4, 4))
            violations = walls + spin(0)
            expected = violations <= 1
            actual = min(energy_change(spin, s, line_incident)[0] for s in flips) >= 0
            assert actual == expected, (left, right, bits, violations)
            if actual:
                boundary += 1
                sectors.add((left, right))
            tested += 1
    assert boundary == 10 and len(sectors) == 4
    # This is a zero-energy movable right-hand wall, as in Prop. 6.23.
    s = frozenset(range(1, 4))
    assert energy_change(lambda i: int(i > 0), s, line_incident)[0] == 0
    print(f'FIDELITY Prop.6.23: configurations={tested}, flips_per_configuration={len(flips)}, boundary={boundary}, asymptotic_sectors={len(sectors)}')


def comb_incident(v):
    i, j = v
    out = {frozenset(((i - 1, j), v)), frozenset((v, (i + 1, j)))}
    if i == 0:
        out.update((frozenset(((0, j - 1), v)), frozenset((v, (0, j + 1)))))
    if v == (0, 0):
        out.add(frozenset((v,)))
    return out


def path_to(v):
    i, j = v
    path = []
    for t in range(min(0, j), max(0, j)):
        path.append(frozenset(((0, t), (0, t + 1))))
    for t in range(min(0, i), max(0, i)):
        path.append(frozenset(((t, j), (t + 1, j))))
    return path


def model_checks():
    for v in product(range(-4, 5), repeat=2):
        total = frozenset(((0, 0),))
        for edge in path_to(v):
            total = total ^ edge
        assert total == frozenset((v,))
        assert len(comb_incident(v)) <= 5
    # Symplectic nondegeneracy and L=L-perp for every Pauli on three sites.
    n = 3
    diagonal = {(a, 0) for a in range(1 << n)}
    def symp(p, q):
        return ((p[0] & q[1]).bit_count() + (p[1] & q[0]).bit_count()) % 2
    paulis = tuple(product(range(1 << n), repeat=2))
    perp = {p for p in paulis if all(symp(p, q) == 0 for q in diagonal)}
    assert perp == diagonal
    assert all(any(symp(p, q) for q in paulis) for p in paulis if p != (0, 0))
    print('MODEL: 81 exact path-generation identities; max incident generators=5; 3-site Lagrangian/symplectic test=OK')


def comb_boundary_checks():
    window = tuple(product(range(-1, 3), range(-1, 2)))
    flips = tuple(subsets(window))
    tests = 0
    max_affected = 0
    for ks in product(range(3), repeat=3):
        def spin(v):
            i, j = v
            k = ks[j + 1] if -1 <= j <= 1 else 0
            return int(i > k)
        for s in flips:
            change, affected = energy_change(spin, s, comb_incident)
            assert len(affected) <= 5 * len(s)
            max_affected = max(max_affected, len(affected))
            assert change >= 0, (ks, s, change)
            # Independent direct local-Hamiltonian difference over all incident
            # terms, including unaffected terms that must cancel exactly.
            incident = set().union(*(comb_incident(v) for v in s)) if s else set()
            direct = sum(parity(spin(v) ^ (v in s) for v in g) - parity(spin(v) for v in g) for g in incident)
            assert direct == change
            # Count total walls on each affected row, using unchanged tails.
            horizontal = 0
            for j in {v[1] for v in s}:
                k = ks[j + 1]
                lo = min([v[0] for v in s if v[1] == j] + [k]) - 1
                hi = max([v[0] for v in s if v[1] == j] + [k]) + 1
                d = sum((spin((i, j)) ^ ((i, j) in s)) !=
                        (spin((i + 1, j)) ^ ((i + 1, j) in s)) for i in range(lo, hi + 1))
                assert d >= 1 and d % 2 == 1
                horizontal += d - 1
            spine_pin = change - horizontal
            assert spine_pin >= 0
            tests += 1
    assert tests == 27 * 4096
    print(f'COMB BOUNDARY: wall_patterns=27, flips_per_pattern=4096, exact_tests={tests}, min_energy_change=0, max_affected={max_affected}')


def orbit_and_cylinder_checks():
    rows = tuple(range(-2, 3))
    patterns = tuple(product(range(3), repeat=len(rows)))
    checked = 0
    for k in patterns:
        for ell in patterns:
            km, lm = dict(zip(rows, k)), dict(zip(rows, ell))
            s = frozenset((i, j) for j in rows for i in range(min(km[j], lm[j]) + 1, max(km[j], lm[j]) + 1))
            assert len(s) == sum(abs(a - b) for a, b in zip(k, ell))
            def x(v): return int(v[0] > km.get(v[1], 0))
            def y(v): return int(v[0] > lm.get(v[1], 0))
            assert energy_change(x, s, comb_incident)[0] == 0
            assert all((x(v) ^ (v in s)) == y(v) for v in product(range(-1, 4), range(-3, 4)))
            checked += 1
    # Basic cylinders at the all-0 and all-1 wall sequences, all subsets of
    # this finite observation window. Modify exactly the finitely observed rows.
    sites = tuple(product(range(-1, 3), range(-1, 2)))
    cylinders = 0
    for observed in subsets(sites):
        s = frozenset((1, j) for _, j in observed)
        x = lambda v: int(v[0] > 0)
        y = lambda v: int(v[0] > 1)
        assert all((x(v) ^ (v in s)) == y(v) for v in observed)
        assert all((y(v) ^ (v in s)) == x(v) for v in observed)
        assert energy_change(x, s, comb_incident)[0] == 0
        assert energy_change(y, s, comb_incident)[0] == 0
        cylinders += 1
    print(f'ORBITS: exact wall-shift pairs={checked}, zero_energy_moves=OK; mutual cylinder approximations={cylinders}')
    print('INFINITE SCOPE: finite checks do not prove the infinite boundary or non-T0 assertion; those require the stated general argument and future Lean probe.')


def main():
    fidelity_biased_ising()
    model_checks()
    comb_boundary_checks()
    orbit_and_cylinder_checks()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

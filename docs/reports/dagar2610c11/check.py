#!/usr/bin/env python3
"""Exact Mackey double-coset check for Dagar, arXiv:2610.15762, Conjecture 1.1.

All matrices are over Z/(p^2). Character values are encoded by exponents of a
primitive (p*(p-1))-th root of unity. No floating point or symbolic rounding.
B-orbits on G/B are double cosets; Schreier generators give their stabilizers.
The code is independent of the proposed transvection/commutator proof.
"""
from collections import Counter, deque
from itertools import product
from math import prod


def identity(n):
    return tuple(int(i == j) for i in range(n) for j in range(n))


def mul(a, b, n, q):
    return tuple(sum(a[i*n+k]*b[k*n+j] for k in range(n)) % q
                 for i in range(n) for j in range(n))


def inverse(a, n, q, p):
    rows = [list(a[i*n:(i+1)*n])+list(identity(n)[i*n:(i+1)*n]) for i in range(n)]
    for j in range(n):
        k = next(k for k in range(j, n) if rows[k][j] % p)
        rows[k], rows[j] = rows[j], rows[k]
        scale = pow(rows[j][j], -1, q)
        rows[j] = [x*scale % q for x in rows[j]]
        for i in range(n):
            if i != j:
                c = rows[i][j]
                rows[i] = [(x-c*y) % q for x, y in zip(rows[i], rows[j])]
    ans = tuple(x for row in rows for x in row[n:])
    assert mul(a, ans, n, q) == identity(n)
    return ans


def canonical(a, n, q, p):
    """Canonical right B-coset by upper triangular column operations."""
    a = list(a)
    used = []
    for j in range(n):
        for k, pivot in enumerate(used):
            c = a[pivot*n+j]
            for i in range(n):
                a[i*n+j] = (a[i*n+j]-c*a[i*n+k]) % q
        pivot = next(i for i in range(n) if i not in used and a[i*n+j] % p)
        scale = pow(a[pivot*n+j], -1, q)
        for i in range(n):
            a[i*n+j] = a[i*n+j]*scale % q
        used.append(pivot)
    return tuple(a)


def generators(n, q, primitive, upper_only=False):
    ans = []
    for i in range(n):
        a = list(identity(n))
        a[i*n+i] = primitive
        ans.append(tuple(a))
    for i in range(n):
        for j in range(n):
            if i != j and (not upper_only or i < j):
                a = list(identity(n))
                a[i*n+j] = 1
                ans.append(tuple(a))
    return ans


def is_borel(a, n, p):
    return (all(a[i*n+j] == 0 for i in range(n) for j in range(i))
            and all(a[i*n+i] % p for i in range(n)))


def solve(p, n):
    q = p*p
    m = p*(p-1)
    primitive = 3 if p == 2 else 2
    logs = {pow(primitive, k, q): k for k in range(m)}
    assert len(logs) == m
    eye = identity(n)
    gg = generators(n, q, primitive)
    bg = generators(n, q, primitive, True)
    flags, queue = {eye}, deque([eye])
    while queue:
        a = queue.popleft()
        for s in gg:
            b = canonical(mul(s, a, n, q), n, q, p)
            if b not in flags:
                flags.add(b)
                queue.append(b)
    group_order = p**(n*n)*prod(p**n-p**i for i in range(n))
    borel_order = m**n*q**(n*(n-1)//2)
    assert len(flags)*borel_order == group_order
    if n == 3:
        # Source, Section 3.1, the displayed dimension of B(chi1,chi2,chi3).
        assert len(flags) == p**3*(p+1)*(p*p+p+1)
    unseen = set(flags)
    data = []
    while unseen:
        root = min(unseen)
        root_inv = inverse(root, n, q, p)
        trans = {root: eye}
        queue = deque([root])
        constraints = set()
        while queue:
            a = queue.popleft()
            for s in bg:
                b = canonical(mul(s, a, n, q), n, q, p)
                sb = mul(s, trans[a], n, q)
                if b not in trans:
                    trans[b] = sb
                    queue.append(b)
                # Schreier stabilizer generator of root B, expressed in B.
                h = mul(inverse(trans[b], n, q, p), sb, n, q)
                k = mul(mul(root_inv, h, n, q), root, n, q)
                assert is_borel(h, n, p) and is_borel(k, n, p)
                constraints.add(tuple((logs[k[i*n+i]]-logs[h[i*n+i]]) % m
                                      for i in range(n)))
        assert borel_order % len(trans) == 0
        unseen.difference_update(trans)
        data.append((root, len(trans), constraints))
    assert sum(size for _, size, _ in data) == len(flags)
    tuples = list(product(range(m), repeat=n))
    norms = {cs: sum(all(sum(c*d for c, d in zip(cs, row)) % m == 0
                           for row in constraints) for _, _, constraints in data)
             for cs in tuples}
    # Independent direct enumeration of B and each B intersection in small groups.
    if n == 2 or p == 2:
        positions = [(i, j) for i in range(n) for j in range(i+1, n)]
        direct = [set() for _ in data]
        roots_inv = [inverse(root, n, q, p) for root, _, _ in data]
        for diag in product(sorted(logs), repeat=n):
            for upper in product(range(q), repeat=len(positions)):
                h = [0]*(n*n)
                for i, x in enumerate(diag):
                    h[i*n+i] = x
                for (i, j), x in zip(positions, upper):
                    h[i*n+j] = x
                h = tuple(h)
                for idx, (root, _, _) in enumerate(data):
                    k = mul(mul(roots_inv[idx], h, n, q), root, n, q)
                    if is_borel(k, n, p):
                        direct[idx].add(tuple((logs[k[i*n+i]]-logs[h[i*n+i]]) % m
                                              for i in range(n)))
        for cs in tuples:
            direct_norm = sum(all(sum(c*d for c, d in zip(cs, row)) % m == 0
                                  for row in constraints) for constraints in direct)
            assert norms[cs] == direct_norm
    for cs, norm in norms.items():
        distinct = len({c % p for c in cs}) == n
        assert (norm == 1) == distinct, (p, n, cs, norm, distinct)
    print(f'p={p} n={n} ell=2 |G|={group_order} |B|={borel_order} '
          f'flags={len(flags)} double_cosets={len(data)} tuples={len(tuples)}')
    print('norm_histogram=' + repr(sorted(Counter(norms.values()).items())))
    for cs in tuples:
        print(f'  character_exponents={cs} mackey_norm={norms[cs]}')
    print('SOURCE_FIDELITY_OK' if n == 3 else 'RANK_TWO_CHECK_OK')
    print('DIRECT_B_INTERSECTION_CROSSCHECK_OK' if n == 2 or p == 2
          else 'SCHREIER_STABILIZERS_EXACT')
    return len(tuples), sum(norm == 1 for norm in norms.values())


def main():
    checked = irreducible = 0
    for p in (2, 3):
        for n in (2, 3):
            count, irr = solve(p, n)
            checked += count
            irreducible += irr
    assert checked == 264 and irreducible == 74
    print(f'TOTAL character_tuples={checked} irreducible={irreducible}')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

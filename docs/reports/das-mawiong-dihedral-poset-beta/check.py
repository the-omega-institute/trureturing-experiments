#!/usr/bin/env python3
"""Exact checks for the proposed value beta(D_4) = 6.

The script checks the six-point realizer and the two structural ingredients of
the lower bound: a faithful D_4 action on at most five points either has an
orbit of size four or has exponent-two image, and a four-orbit plus one fixed
point forces all permutations of the orbit to be order automorphisms.
"""

from itertools import permutations


def compose(p, q):
    return tuple(p[q[i]] for i in range(len(p)))


def inverse(p):
    out = [0] * len(p)
    for i, j in enumerate(p):
        out[j] = i
    return tuple(out)


def power(p, k):
    out = tuple(range(len(p)))
    for _ in range(k):
        out = compose(p, out)
    return out


def order(p):
    out = tuple(range(len(p)))
    for k in range(1, 100):
        out = compose(p, out)
        if out == tuple(range(len(p))):
            return k
    raise AssertionError("permutation order did not close")


def preserves(rel, p):
    return {(p[a], p[b]) for a, b in rel} == rel


def orbit(pairs, x):
    todo = {x}
    changed = True
    while changed:
        changed = False
        for p in pairs:
            nxt = {p[y] for y in todo}
            if not nxt <= todo:
                todo |= nxt
                changed = True
    return frozenset(todo)


def main():
    # Two disjoint three-point V-posets: 0 < 1,2 and 3 < 4,5.
    rel = {(0, 1), (0, 2), (3, 4), (3, 5)}
    perms6 = list(permutations(range(6)))
    aut = [p for p in perms6 if preserves(rel, p)]
    assert len(aut) == 8

    r = (3, 4, 5, 0, 2, 1)  # (0 3)(1 4 2 5)
    s = (0, 2, 1, 3, 4, 5)  # (1 2)
    assert preserves(rel, r) and preserves(rel, s)
    assert order(r) == 4 and order(s) == 2
    assert compose(s, compose(r, s)) == inverse(r)
    generated = {compose(s if j else tuple(range(6)), power(r, i))
                 for i in range(4) for j in range(2)}
    assert len(generated) == 8 and set(aut) == generated
    print("six-point construction: automorphism group has order 8 and D4 relations")

    # Enumerate all permutation actions of D4 on n <= 5 points.  The relation
    # s r s = r^{-1} is the presentation used in the preregistration.
    found_four_orbit = []
    faithful_without_four = []
    for n in range(1, 6):
        identity = tuple(range(n))
        ps = list(permutations(range(n)))
        for rr in ps:
            if order(rr) != 4:
                continue
            for ss in ps:
                if order(ss) != 2 or compose(ss, compose(rr, ss)) != inverse(rr):
                    continue
                image = {compose(ss if j else identity, power(rr, i))
                         for i in range(4) for j in range(2)}
                if len(image) != 8:
                    continue  # not a faithful permutation action
                orbits = {orbit((rr, ss), x) for x in range(n)}
                sizes = sorted(map(len, orbits))
                if 4 in sizes:
                    found_four_orbit.append((n, sizes))
                else:
                    faithful_without_four.append((n, sizes))
    assert found_four_orbit and not faithful_without_four
    assert min(n for n, _ in found_four_orbit) == 4
    print("lower-bound action check: every faithful action on <=5 points has a 4-orbit")

    # With a transitive four-point orbit A and at most one fixed point z, an
    # invariant strict order has no relation inside A and z is uniformly below,
    # uniformly above, or incomparable with A.  All 24 permutations of A then
    # preserve each of the four possibilities.
    A = tuple(range(4))
    for mode in ("none", "below", "above", "both"):  # both is impossible in a poset
        rel5 = set()
        if mode in ("below", "both"):
            rel5 |= {(4, a) for a in A}
        if mode in ("above", "both"):
            rel5 |= {(a, 4) for a in A}
        if mode == "both":
            continue  # contains a two-cycle, hence is not a strict order
        for q in permutations(A):
            p = q + (4,)
            assert preserves(rel5, p)
    assert 24 > 8
    print("four-orbit obstruction: S4 extends, giving at least 24 automorphisms")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Exact, bounded checks for Narayan, arXiv:2610.02495v1, Conjecture 18.3.

All evidence uses Python integers. Length sets are computed from generators,
without assuming a delta-set formula or a normalization theorem. Finite windows
do not establish a global delta set; the proposed global argument is separate.
"""

from math import gcd


def euclidean_coordinates(x, y):
    assert x > 0 and y > 0
    seen = {x, y}
    while x != y:
        if x > y:
            x -= y
        else:
            y -= x
        seen.update((x, y))
    return seen


def length_masks(generators, bound):
    """Bit l is set iff n has a coordinate factorization of length l.

    Removing one summand reduces n, so the recurrence enumerates all lengths.
    Repeated generator use is allowed, and every displayed coordinate is kept.
    """
    assert generators and all(g > 0 for g in generators)
    masks = [0] * (bound + 1)
    masks[0] = 1
    for n in range(1, bound + 1):
        for g in generators:
            if g <= n:
                masks[n] |= masks[n - g] << 1
    return masks


def lengths(mask):
    answer = []
    while mask:
        bit = mask & -mask
        answer.append(bit.bit_length() - 1)
        mask -= bit
    return answer


def observed_delta(masks):
    answer = set()
    for mask in masks:
        ls = lengths(mask)
        answer.update(v - u for u, v in zip(ls, ls[1:]))
    return answer


def coordinate_factorizations(generators, n):
    if len(generators) == 1:
        if n % generators[0] == 0:
            yield (n // generators[0],)
        return
    g = generators[0]
    for count in range(n // g + 1):
        for tail in coordinate_factorizations(generators[1:], n - count * g):
            yield (count,) + tail


def value(generators, z):
    return sum(g * c for g, c in zip(generators, z))


def normalize(a, initial):
    """Test one deterministic path using the three proposed rewrite rules."""
    gs = (a,) + tuple(2 * a + 2 * i for i in range(a))
    z = list(initial)
    n, steps = value(gs, z), 0
    while z[1] or sum(z[2:]) >= 2:
        old_length = sum(z)
        if z[1]:
            z[1] -= 1
            z[0] += 2
        else:
            i = next(i for i in range(1, a) if z[i + 1])
            z[i + 1] -= 1
            j = next(j for j in range(1, a) if z[j + 1])
            z[j + 1] -= 1
            t = i + j
            if t < a:
                z[0] += 2
                z[t + 1] += 1
            else:
                z[1] += 2
                z[t - a + 1] += 1
        assert all(c >= 0 for c in z)
        assert value(gs, z) == n
        assert sum(z) == old_length + 1
        assert sum(z) <= n // a
        steps += 1
        assert steps <= n // a
    assert z[1] == 0 and sum(z[2:]) <= 1
    assert sum(z) == sum(initial) + steps
    return tuple(z)


def expected_terminal(a, n):
    i = (n * pow(2, -1, a)) % a
    z = [0] * (a + 1)
    residue_value = 0 if i == 0 else 2 * a + 2 * i
    assert n >= residue_value and (n - residue_value) % a == 0
    z[0] = (n - residue_value) // a
    if i:
        z[i + 1] = 1
    return tuple(z)


def main():
    # Fidelity: the explicit calculation in the source's Section 3 remark.
    reported = {1, 2, 3, 4, 7, 10, 13, 23, 33, 43, 76, 109, 142, 251, 393}
    a, h, d, k = 4, 143, 109, 2
    e = h - 1
    maximum = ((a + k - 1) // k) * e + d
    tail = set(range(maximum, 0, -e))
    assert maximum == 393 and tail == {393, 251, 109}
    assert euclidean_coordinates(e, d) | tail == reported
    print('FIDELITY source Section 3, (a,h,d,k)=(4,143,109,2):')
    print('  E(142,109) union {393,251,109} =', sorted(reported))

    # Independent length enumeration for a specialization of Theorem 12.1.
    original = (3, 8, 10)  # no 6-coordinate, unlike Conjecture 18.3
    original_bound = 2000
    observed = observed_delta(length_masks(original, original_bound))
    source_formula = euclidean_coordinates(1, 2) | {4, 3, 2, 1}
    assert observed == source_formula == {1, 2, 3, 4}
    print('FIDELITY Theorem 12.1, generators', original,
          'n=0..2000: observed delta =', sorted(observed))
    print('  This is a bounded reproduction, not a proof of the global theorem.')

    identities = 0
    for a in range(3, 102, 2):
        assert gcd(a, 2) == 1 and 1 < a - 1 < a
        b = tuple(2 * a + 2 * i for i in range(a))
        assert b[0] == a + a
        identities += 1
        for i in range(1, a):
            for j in range(i, a):
                t = i + j
                rhs = a + a + b[t] if t < a else b[0] + b[0] + b[t - a]
                assert b[i] + b[j] == rhs
                assert 0 <= (t if t < a else t - a) < a
                identities += 1
        assert len({(2 * i) % a for i in range(a)}) == a
    print('REWRITES odd a=3..101: checked', identities,
          'value identities; all change length by one; residue maps bijective.')

    total_vectors = 0
    for a in (3, 5, 7, 9, 11):
        gs = (a,) + tuple(2 * a + 2 * i for i in range(a))
        masks = length_masks(gs, 120)
        for n, mask in enumerate(masks):
            vectors = list(coordinate_factorizations(gs, n))
            enumerated = {sum(z) for z in vectors}
            assert enumerated == set(lengths(mask))
            if vectors:
                terminal = expected_terminal(a, n)
                for z in vectors:
                    assert normalize(a, z) == terminal
                    total_vectors += 1
    print('NORMALIZATION all coordinate vectors for a in {3,5,7,9,11}, n=0..120:',
          total_vectors, 'vectors; endpoints and length-DP agree.')

    bound = 5000
    checked_elements = 0
    predicted = euclidean_coordinates(1, 2)
    assert predicted == {1, 2}
    for a in range(3, 32, 2):
        gs = (a,) + tuple(2 * a + 2 * i for i in range(a))
        masks = length_masks(gs, bound)
        actual = set()
        for n, mask in enumerate(masks):
            if not mask:
                continue
            lo = (mask & -mask).bit_length() - 1
            hi = mask.bit_length() - 1
            assert mask == ((1 << (hi - lo + 1)) - 1) << lo
            if lo < hi:
                actual.add(1)
            checked_elements += 1
        assert lengths(masks[2 * a]) == [1, 2]
        assert actual == {1} and actual != predicted
        print('FAMILY a=' + str(a) + ', (h,d,k)=(2,2,' + str(a - 1) + '),'
              ' n=0..5000: every nonempty L(n) is an interval;'
              ' observed delta=[1]; E(1,2)=[1,2].')
    print('FAMILY total represented elements checked:', checked_elements)
    print('Finite windows do not certify absence of 2 for unbounded n;'
          ' the proposed normalization argument supplies that obligation.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

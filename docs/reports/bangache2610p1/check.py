#!/usr/bin/env python3
"""Exact checks for Bangachev, arXiv:2610.27147v1, Open Problem 1.

The obstruction concerns exact exponent length, not an upper bound on length.
All evidence uses integers and Fraction; no numerical limit approximation is used.
"""

from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, product
from math import comb, prod


def parity(bits):
    return sum(bits) % 2


def bernoulli_mass(bits, eta):
    return prod((eta if b else 1 - eta) for b in bits)


def xor_probability(length, eta):
    return sum((bernoulli_mass(bits, eta)
                for bits in product(range(2), repeat=length) if parity(bits)), Q(0))


def collision_free_probability(n, m):
    return Q(prod(n - i for i in range(m)), n**m) if m <= n else Q(0)


def fidelity():
    cases = 0
    for length in range(1, 9):
        for eta in [Q(0), Q(1, 16), Q(1, 5), Q(1, 2), Q(1)]:
            p = xor_probability(length, eta)
            assert p == (1 - (1 - 2 * eta)**length) / 2
            assert p <= length * eta  # the union bound used by Lemma 1
            cases += 1
    cancellations = 0
    # Reproduce the source's cycle-to-test cancellation for every tiny matrix,
    # including both sparsities in dimension 2, and every nonempty cycle.
    for n in range(1, 4):
        for k in range(1, n + 1):
            rows = [sum(1 << j for j in support)
                    for support in combinations(range(n), k)]
            for m in range(1, 4):
                for a in product(rows, repeat=m):
                    for mask in range(1, 1 << m):
                        indices = [i for i in range(m) if mask >> i & 1]
                        row_sum = 0
                        for i in indices:
                            row_sum ^= a[i]
                        if row_sum:
                            continue
                        for secret in range(1 << n):
                            signal = [parity([(a[i] & secret).bit_count()])
                                      for i in range(m)]
                            for errors in product(range(2), repeat=m):
                                y = [(signal[i] + errors[i]) % 2 for i in range(m)]
                                assert parity([y[i] for i in indices]) == parity(
                                    [errors[i] for i in indices])
                                cancellations += 1
    print(f"FIDELITY: Lemma 1 XOR law/union bound: {cases} exact cases")
    print(f"FIDELITY: cycle signal cancellation: {cancellations} exact cases")


def occupancy_and_cycles():
    matrices = subsets = 0
    for n in range(1, 5):
        for m in range(1, 7):
            collision_sum = 0
            no_collision = 0
            for a in product(range(n), repeat=m):
                c = Counter(a)
                d = sum(comb(count, 2) for count in c.values())
                collision_sum += d
                no_collision += int(d == 0)
                maximum = 2 * sum(count // 2 for count in c.values())
                assert maximum <= 2 * d
                actual_maximum = 0
                for mask in range(1 << m):
                    selected = Counter(a[i] for i in range(m) if mask >> i & 1)
                    zero_sum = all(count % 2 == 0 for count in selected.values())
                    if zero_sum:
                        length = mask.bit_count()
                        actual_maximum = max(actual_maximum, length)
                        assert length <= 2 * d
                    subsets += 1
                assert maximum == actual_maximum
                matrices += 1
            assert Q(collision_sum, n**m) == Q(comb(m, 2), n)
            assert Q(no_collision, n**m) == collision_free_probability(n, m)
    print(f"OCCUPANCY: {matrices} ordered matrices; {subsets} row subsets")
    print("OCCUPANCY: exact birthday probability, E[D]=choose(m,2)/n, max L=2 sum floor(C/2)<=2D")


def tester():
    cases = 0
    # For each A, enumerate every secret and error vector and every H0 response.
    # A deterministic fallback outputs 0 if no duplicate exists.
    for n, m in [(2, 3), (3, 4), (4, 3)]:
        for eta in [Q(1, 25), Q(1, 5), Q(1, 2)]:
            h0 = Q(0)
            h1 = Q(0)
            for a in product(range(n), repeat=m):
                pairs = [(i, j) for i in range(m) for j in range(i + 1, m)
                         if a[i] == a[j]]
                pair = pairs[0] if pairs else None
                def outputs_zero(y):
                    return pair is None or parity([y[pair[0]], y[pair[1]]]) == 1
                h0 += Q(sum(outputs_zero(y) for y in product(range(2), repeat=m)),
                        n**m * 2**m)
                for secret in product(range(2), repeat=n):
                    for errors in product(range(2), repeat=m):
                        y = [(secret[a[i]] + errors[i]) % 2 for i in range(m)]
                        if outputs_zero(y):
                            h1 += bernoulli_mass(errors, eta) / (n**m * 2**n)
            p_e = 1 - collision_free_probability(n, m)
            q = 2 * eta * (1 - eta)
            assert h0 == 1 - p_e / 2
            assert h1 == 1 - p_e + p_e * q
            assert h0 - h1 == p_e * (Q(1, 2) - q)
            cases += 1
    print(f"TESTER: {cases} full exact H0/H1 enumerations; advantage=P(E)(1/2-2 eta(1-eta))")


def exact_exponent_family():
    # n=t^5, m=t^3, eta=t^-2: no rounding and no irrational powers.
    # The elementary birthday bound prod(1-i/n)<=1/(1+sum(i/n))
    # gives P(E)>=a/(1+a), a=choose(m,2)/n.
    for t in [2, 3, 4, 10, 100, 10000]:
        n, m, eta = t**5, t**3, Q(1, t**2)
        a = Q(comb(m, 2), n)
        assert a == Q(t**3 - 1, 2 * t**2)
        assert a >= Q(t, 4)
        duplicate_lower = a / (1 + a)
        advantage_lower = duplicate_lower * (Q(1, 2) - 2 * eta * (1 - eta))
        assert advantage_lower >= Q(1, 24) > Q(1, 100)
        # Any L >= n^(3/10)=t^(3/2) requires 2D >= t^(3/2).
        # Markov <= 2 E[D]/t^(3/2) <= t^(-1/2).
        # Square positive quantities to check these irrational comparisons exactly.
        assert (2 * a)**2 <= t**2
        assert t**3 < t**4  # (t^(3/2))^2 < (n^(2/5))^2
        if t <= 4:
            p0 = collision_free_probability(n, m)
            assert p0 <= 1 / (1 + a)
        print(f"FAMILY: t={t}; n={n}; m={m}; eta={eta}; E[D]={a}; advantage_lower={advantage_lower}")

    families = 0
    for denominator in range(3, 21):
        for numerator in range(1, denominator):
            beta = Q(numerator, denominator)
            if not Q(1, 2) < beta < 1:
                continue
            # Choose a rational delta strictly between 2 beta-1 and beta.
            delta = (3 * beta - 1) / 2
            alpha = (2 * beta - 1 + delta) / 2
            assert 0 < 2 * beta - 1 < alpha < delta < beta < 1
            # For n=t^denominator, m=t^numerator, E[D] is this exact rational.
            for t in range(2, 8):
                n, m = t**denominator, t**numerator
                assert m**denominator == n**numerator
                expected = Q(comb(m, 2), n)
                assert 2 * expected <= t**(2 * numerator - denominator)
            families += 1
    print(f"FAMILY: {families} rational exponent pairs; strict gaps 2 beta-1<alpha<delta<beta")
    print("SCOPE: exact exponent n^(delta+o(1)) only; a 2-cycle satisfies an upper-bound interpretation")


def main():
    fidelity()
    occupancy_and_cycles()
    tester()
    exact_exponent_family()
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()

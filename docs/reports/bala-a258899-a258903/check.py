#!/usr/bin/env python3
"""Exact finite checks for Bala's OEIS A258899 and A258903 conjectures.

Standard library only. These checks are evidence on the stated finite domains;
they do not prove the unbounded conjectures. Sources and interpretation: README.md.
"""
from fractions import Fraction
from math import comb, factorial, prod

MAX_INDEX = 240
MAX_MODULUS = 60
FACT = [factorial(n) for n in range(MAX_INDEX + 1)]
PREFIXES = {
    (2, 1): [1, 2, 6, 10, 42, 34, 786, 130, 17058, 81154, 545346,
             2050, 102457218, 8194, 1141636866, 72648608770,
             648648065538, 131074, 111258180895746, 524290],
    (2, 3): [1, 6, 30, 78, 426, 582, 12450, 4758, 407010, 2218182,
             19172370, 360438, 4755166050, 3213222, 85631151090,
             5099958831318, 54483404779650],
}


def source_sequence(u, v, bound):
    """The source finite sum, using exact rational division before integrality."""
    a = [1]
    for n in range(1, bound + 1):
        total = sum((Fraction(u**d * v**(n // d) * FACT[n],
                              FACT[d] * FACT[n // d])
                     for d in range(1, n + 1) if n % d == 0), Fraction(0))
        assert total.denominator == 1, (u, v, n, total)
        a.append(total.numerator)
    return a


def egf_coefficients(u, v, bound):
    """Locally finite expansion of 1 + sum_d u^d/d! (exp(v*x^d)-1)."""
    coeff = [Fraction(1)] + [Fraction(0) for _ in range(bound)]
    for d in range(1, bound + 1):
        for e in range(1, bound // d + 1):
            coeff[d * e] += Fraction(u**d * v**e, FACT[d] * FACT[e])
    return coeff


def phi(m):
    from math import gcd
    return sum(gcd(j, m) == 1 for j in range(1, m + 1))


def inverse_stirling(a):
    rows = [[1]]
    weights = []
    for n in range(len(a)):
        if n:
            prev = rows[-1]
            rows.append([0] + [
                (prev[r - 1] if r - 1 < len(prev) else 0)
                + r * (prev[r] if r < len(prev) else 0)
                for r in range(1, n + 1)])
        weights.append(a[n] - sum(weights[r] * rows[n][r] for r in range(n)))
    return weights


def cycle_weights(a):
    """Unique w_n from F=exp(sum_n w_n*x^n/n), using F'=F*(log F)'."""
    f = [Fraction(a[n], FACT[n]) for n in range(len(a))]
    w = [Fraction(0)]
    for n in range(1, len(a)):
        w.append(n * f[n] - sum(w[k] * f[n - k] for k in range(1, n)))
    return w


def main():
    sources = {uv: source_sequence(*uv, MAX_INDEX) for uv in PREFIXES}
    for uv, prefix in PREFIXES.items():
        assert sources[uv][:len(prefix)] == prefix
        coeff = egf_coefficients(*uv, 40)
        assert [coeff[n] * FACT[n] for n in range(41)] == sources[uv][:41]
        print(f"FIDELITY u={uv[0]} v={uv[1]}: OEIS prefix n=0..{len(prefix)-1}; "
              "independent e.g.f. expansion n=0..40: OK")
    a, b = sources[(2, 1)], sources[(2, 3)]
    mod11 = [1, 2, 6, 10, 9, 1, 5, 9, 8, 7, 10, 4, 6, 10, 7, 1, 0,
             9, 5, 8, 3, 4, 6, 10, 7, 1, 0, 9, 5, 8, 3, 4, 6, 10, 7, 1,
             0, 9, 5, 8, 3]
    mod7 = [6, 2, 1, 6, 1, 4, 5, 2, 1, 0, 1, 5, 5, 2, 1, 0, 1, 5, 5,
            2, 1, 0, 1, 5, 5, 2, 1, 0, 1, 5]
    assert [x % 11 for x in a[:len(mod11)]] == mod11
    # A258903's comment starts at n=1; the entry's offset is zero.
    assert [x % 7 for x in b[1:len(mod7)+1]] == mod7
    assert all((a[n+10]-a[n]) % 11 == 0 for n in range(11, MAX_INDEX-9))
    assert all((b[n+6]-b[n]) % 7 == 0 for n in range(7, MAX_INDEX-5))
    print("FIDELITY source residue lists: A258899 starts n=0; A258903 starts n=1: OK")
    print("FIDELITY displayed periods: A258899 mod 11 from n=11; A258903 mod 7 from n=7, through n+T=240: OK")

    proper_count = 0
    for n in range(1, MAX_INDEX+1):
        for d in range(1, n+1):
            if n % d:
                continue
            e = n // d
            q = Fraction(FACT[n], FACT[d] * FACT[e])
            assert q.denominator == 1
            if d >= 2 and e >= 2:
                assert d+e <= n and 2*(n-d-e) >= n-4
                assert q == comb(n, d) * prod(range(e+1, n-d+1))
                for m in range(1, MAX_MODULUS+1):
                    if n >= 2*m+4:
                        assert n-d-e >= m
                        assert q.numerator % m == 0
                        proper_count += 1
    print(f"FACTORIAL proper-term quotient/product and divisibility checks: {proper_count}: OK")

    tail_checks = period_checks = power_checks = 0
    for u in range(1, 5):
        for v in range(1, 5):
            seq = sources.get((u, v))
            if seq is None:
                seq = source_sequence(u, v, MAX_INDEX)
            for m in range(1, MAX_MODULUS+1):
                t = phi(m)
                assert t > 0
                for n in range(2*m+4, MAX_INDEX+1):
                    assert (seq[n] - (u*v**n + v*u**n)) % m == 0
                    tail_checks += 1
                    if n+t <= MAX_INDEX:
                        assert (seq[n+t] - seq[n]) % m == 0
                        period_checks += 1
    for z in range(21):
        for m in range(1, MAX_MODULUS+1):
            t = phi(m)
            for n in range(m, MAX_INDEX-t+1):
                assert (pow(z, n+t, m)-pow(z, n, m)) % m == 0
                power_checks += 1
    print(f"TAIL u,v=1..4, m=1..{MAX_MODULUS}, 2m+4<=n<=240: {tail_checks}: OK")
    print(f"TOTIENT SHIFT same parameters, n+phi(m)<=240: {period_checks}: OK")
    print(f"POWER z=0..20, m=1..{MAX_MODULUS}, m<=n, n+phi(m)<=240: {power_checks}: OK")

    da, db = inverse_stirling(a[:17]), inverse_stirling(b[:17])
    assert Fraction(da[3], FACT[3]) == Fraction(-2, 3)
    assert Fraction(db[4], FACT[4]) == Fraction(21, 2)
    assert da[8] == 421176 and da[8] % 5 == 1
    assert db[8] == 5363064 and db[8] % 5 == 4
    wa, wb = cycle_weights(a[:7]), cycle_weights(b[:7])
    assert wa[4] == Fraction(23, 3)
    assert wb[5] == Fraction(-3269, 4)
    print("PRIOR-ART applicability: inverse-Stirling g3(A258899)=-2/3; g4(A258903)=21/2: OK")
    print("DANESH Theorem 3.1(ii) fails at p=5,r=8: d8=421176 and 5363064; valuation 0 < ell_5(8)-1=1: OK")
    print("KALLAT integral-cycle condition fails: w4(A258899)=23/3; w5(A258903)=-3269/4: OK")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()

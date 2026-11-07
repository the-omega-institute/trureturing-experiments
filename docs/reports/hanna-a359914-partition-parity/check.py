#!/usr/bin/env python3
"""Exact finite checks of Hanna's first COMMENTS conjecture in OEIS A359914.

Source: Paul D. Hanna, https://oeis.org/A359914, revision #43.
The source's literal bilateral equation is evaluated after distributing its
two monomials, cancelling negative x powers before entering Z[[x]].
No claim of an all-index proof is made by this finite computation.
"""

DEGREE = 40
SOURCE_INITIAL = (
    1, 2, 4, 30, 154, 1078, 7046, 50766, 364268, 2713444,
    20384884, 155954760, 1204192106, 9400024042, 73945396990,
    586088682472, 4673927031694, 37484566094970, 302098932029282,
    2445538771089012, 19875632898821430, 162118004651048048,
    1326658157736876148,
)


# Paul D. Hanna, OEIS b359914.txt, n=0..40.
BFILE_0_TO_40 = (1, 2, 4, 30, 154, 1078, 7046, 50766, 364268, 2713444, 20384884, 155954760, 1204192106, 9400024042, 73945396990, 586088682472, 4673927031694, 37484566094970, 302098932029282, 2445538771089012, 19875632898821430, 162118004651048048, 1326658157736876148, 10888864938350905310, 89617474864220375854, 739425854465095436588, 6115105164625023763944, 50681086929018887228996, 420875603396783831306172, 3501608641308049595422610, 29183246713857123074998088, 243614499095924697626586516, 2036723690359767016733323642, 17052203515644760047745909694, 142959098195778714064377617964, 1200031429547979481724923888438, 10085382700043182463316054591634, 84856083659620945250788583747266, 714724295766787397279090767014686, 6026084006201926816042662074111962, 50857007480140862508314997400783156)

def mul(a, b, degree):
    return [sum(a[j] * b[n-j] for j in range(n+1))
            for n in range(degree+1)]


def inverse(a, degree):
    assert a[0] == 1
    b = [1] + [0] * degree
    for n in range(1, degree+1):
        b[n] = -sum(a[j] * b[n-j] for j in range(1, n+1))
    assert mul(a, b, degree) == [1] + [0] * degree
    return b


def power(a, exponent, degree):
    """Integer powers in the truncated unit group; integers throughout."""
    base = list(a[:degree+1])
    if exponent < 0:
        base, exponent = inverse(base, degree), -exponent
    result = [1] + [0] * degree
    while exponent:
        if exponent % 2:
            result = mul(result, base, degree)
        base = mul(base, base, degree)
        exponent //= 2
    return result


def shifts(k):
    lower = k * (3*k-1) // 2
    upper = k * (3*k+1) // 2
    assert 2*lower == k*(3*k-1)
    assert upper == lower+k
    assert min(lower, upper) >= abs(k)
    return lower, upper


def rhs_original(a, degree, window=None):
    """Sum (-1)^k x^{k(3k-1)/2} A^{3k}(1+x^k A).

    This is the same finite coefficient window and the same distributed
    summand as conventions.lean, with no reindexing or parity shortcuts.
    """
    window = degree if window is None else window
    out = [0] * (degree+1)
    for k in range(-window, window+1):
        lo, hi = shifts(k)
        sign = -1 if k % 2 else 1
        for shift, exponent in [(lo, 3*k), (hi, 3*k+1)]:
            if shift <= degree:
                term = power(a, exponent, degree-shift)
                for j, value in enumerate(term):
                    out[shift+j] += sign*value
    return out


def solve(degree):
    a = [1]
    for d in range(1, degree+1):
        a.append(0)
        # The k=0 term is 1+A; all other shifts are positive.
        a[d] = -rhs_original(a, d)[d]
        assert rhs_original(a, d) == [2] + [0]*d
    return a


def rhs_paired(a, degree):
    out = [0] * (degree+1)
    for k in range(-degree, degree+1):
        shift = shifts(k)[0]
        if shift <= degree:
            sign = -1 if k % 2 else 1
            u = power(a, 3*k, degree-shift)
            v = power(a, 1-3*k, degree-shift)
            for j in range(degree-shift+1):
                out[shift+j] += sign*(u[j]+v[j])
    return out


def partitions(degree):
    # Independent unrestricted-partition count by part multiplicities.
    p = [1] + [0]*degree
    for part in range(1, degree+1):
        for n in range(part, degree+1):
            p[n] += p[n-part]
    return p


def euler_product(degree):
    out = [1] + [0]*degree
    for part in range(1, degree+1):
        factor = [0]*(degree+1)
        factor[0], factor[part] = 1, -1
        out = mul(out, factor, degree)
    return out


def watson_product(a, degree):
    # Independent source FORMULA (2), truncated exactly in Z[[x]].
    out = [1] + [0]*degree
    for n in range(1, degree+2):
        for shift, exponent, sign in [
            (n, 0, -1), (n-1, 1, 1), (n, -1, 1),
            (2*n-1, 2, -1), (2*n-1, -2, -1),
        ]:
            if shift > degree:
                continue
            factor = [1] + [0]*degree
            for j, value in enumerate(power(a, exponent, degree-shift)):
                factor[shift+j] += sign*value
            out = mul(out, factor, degree)
    return out


def main():
    a = solve(DEGREE)
    assert tuple(a[:len(SOURCE_INITIAL)]) == SOURCE_INITIAL
    assert tuple(a) == BFILE_0_TO_40
    # Different windows are equal, including each coefficient's own window.
    rhs = rhs_original(a, DEGREE)
    assert rhs == [2] + [0]*DEGREE
    assert rhs_original(a, DEGREE, DEGREE+7) == rhs
    assert all(rhs_original(a, n)[n] == rhs[n] for n in range(DEGREE+1))
    assert rhs_paired(a, DEGREE) == rhs
    assert watson_product(a, DEGREE) == rhs
    p = partitions(DEGREE)
    assert p[:10] == [1, 1, 2, 3, 5, 7, 11, 15, 22, 30]
    e = [0]*(DEGREE+1)
    for k in range(-DEGREE, DEGREE+1):
        shift = shifts(k)[0]
        if shift <= DEGREE:
            e[shift] += -1 if k % 2 else 1
    assert e == euler_product(DEGREE)
    assert mul(e, p, DEGREE) == [1] + [0]*DEGREE
    assert all(a[n] % 2 == 0 and (a[n]//2-p[n]) % 2 == 0
               for n in range(1, DEGREE+1))
    b = [0] + [value//2 for value in a[1:]]
    one_plus_b = [1] + b[1:]
    expected = [2*value for value in mul(one_plus_b, e, DEGREE)]
    assert all((u-v) % 4 == 0 for u, v in zip(rhs, expected))
    # Test the proposed pairing on other units, independently of the solution.
    for test_degree in [1, 2, 8, 16]:
        for seed in range(5):
            u = [1] + [2*((seed+3)*j % 9-4) for j in range(1, test_degree+1)]
            h = [1] + [value//2 for value in u[1:]]
            target = mul(h, e[:test_degree+1], test_degree)
            original = rhs_original(u, test_degree)
            assert rhs_paired(u, test_degree) == original
            assert all((v-2*w) % 4 == 0 for v, w in zip(original, target))
    print('Arithmetic: Python integers only; no floating-point operations')
    print(f'Fidelity: source NAME/FORMULA (1), EXAMPLE and all 23 listed terms reproduced')
    print(f'Fidelity: source b-file coefficients n=0..{DEGREE} reproduced exactly')
    print('a(0..9) =', a[:10])
    print('p(0..9) =', p[:10])
    print(f'Original bilateral equation and enlarged/coefficientwise windows: n=0..{DEGREE} OK')
    print(f'Reindexed pairing and source Watson product FORMULA (2): n=0..{DEGREE} OK')
    print(f'Euler product and ordinary-partition reciprocal: n=0..{DEGREE} OK')
    print(f'Evenness and a(n)/2 == p(n) modulo 2: n=1..{DEGREE} OK')
    print('Equation-specific modulo-four pairing: solution and 20 additional units OK')
    print('Finite checks support the reading and outline; no all-index proof or Lean probe')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

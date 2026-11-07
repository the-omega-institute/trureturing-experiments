#!/usr/bin/env python3
"""Exact finite checks for Han--Pedon, arXiv:2610.36132v1, Conjecture C.

No floating point, external input, repository imports, or Lean proof is used.
The continuant data are transcribed from the n+1-shift H-fraction package
of Han--Pedon 2502.05993v2 (the separate golden cycle is included).
These checks do not prove any unbounded quantified assertion.
"""
from math import comb


def trim(a):
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return a


def add(a, b):
    return trim([(a[i] if i < len(a) else 0) +
                 (b[i] if i < len(b) else 0) for i in range(max(len(a), len(b)))])


def scale(a, s):
    return trim([s * x for x in a])


def sub(a, b):
    return add(a, scale(b, -1))


def mul(a, b, cap=None):
    size = len(a) + len(b) - 1 if a and b else 0
    if cap is not None:
        size = min(size, cap)
    out = [0] * size
    for i, x in enumerate(a[:size]):
        if x:
            for j, y in enumerate(b[:size-i]):
                out[i+j] += x*y
    return trim(out)


def monomial(k, c=1):
    return [0]*k + [c] if c else []


def coeff(a, k):
    return a[k] if 0 <= k < len(a) else 0


def inv(a, cap):
    assert coeff(a, 0) == 1
    b = [1]
    for k in range(1, cap):
        b.append(-sum(coeff(a, i)*b[k-i] for i in range(1, k+1)))
    return trim(b)


def power(a, k, cap):
    out = [1]
    while k:
        if k % 2:
            out = mul(out, a, cap)
        a = mul(a, a, cap)
        k //= 2
    return out


IDENTITY = [[[1], []], [[], [1]]]


def matmul(a, b, cap=None):
    return [[add(mul(a[i][0], b[0][j], cap),
                 mul(a[i][1], b[1][j], cap)) for j in range(2)] for i in range(2)]


def matadd(a, b):
    return [[add(a[i][j], b[i][j]) for j in range(2)] for i in range(2)]


def matscale(a, p, cap):
    return [[mul(a[i][j], p, cap) for j in range(2)] for i in range(2)]


def matpower(a, k, cap):
    out = IDENTITY
    while k:
        if k % 2:
            out = matmul(out, a, cap)
        a = matmul(a, a, cap)
        k //= 2
    return out


def truncate(a, cap):
    return [[trim(x[:cap]) for x in row] for row in a]


def cycle(n):
    if n == 1:
        return [(0, -1, [1, 1]), (0, 1, [1, 1]), (1, -1, [1, 1, -1])]
    t = add(mul(monomial(1), [1]*n),
            mul(add([1], monomial(n)), [1, -1]))
    states = [('v2', 0), ('v3', 0)]
    for i in range(n-2):
        states += [('w0', i), ('w1', i), ('w2', i)]
    states += [('w0', n-2)]
    for i in range(n-1):
        states += [('u0', i), ('u1', i), ('u2', i)]
    states += [('v0', 0), ('v1', 0)]
    out = []
    for s, i in states:
        if s in ('v0', 'v2'):
            out.append((n-1, 1 if s == 'v0' else -1, add(t, monomial(n+1))))
        elif s == 'v1':
            out.append((n, -1, t))
        elif s == 'v3':
            out.append((0, 1, [1]))
        elif s in ('u0', 'w0'):
            out.append((n-i-2, -1, sub([1]*(n-i), monomial(1))
                        if s == 'w0' else [1]*(n-i)))
        elif s in ('u1', 'w1'):
            out.append((i, -1, sub([1]*(i+2), monomial(1))
                        if s == 'u1' else [1]*(i+2)))
        else:
            out.append((0, -1, [1] if s == 'u2' and i == n-2 else [1, -1]))
    return out


def period(n):
    return (2 if n % 2 == 0 else 4)*n*(n+1)


def blocks(n):
    data = cycle(n)*(1 if n % 2 == 0 else 2)
    u, r = IDENTITY, IDENTITY
    for p, (k, v, d) in enumerate(data):
        previous_k = data[p-1][0]
        u = matmul([[d, monomial(previous_k+k+2, -v)], [[1], []]], u)
        reverse_d = [coeff(d, k+1-i) for i in range(k+2)]
        r = matmul([[trim(reverse_d), [-v]], [[1], []]], r)
    return data, u, r


def moments(n, length):
    b = sub(mul(add([1], monomial(n)), [1, -1]), mul(monomial(1), [1]*n))
    f = [1]
    for k in range(1, length):
        f.append(-sum(f[i]*f[k-1-i] for i in range(k)) -
                 sum(coeff(b, i)*f[k-i] for i in range(1, k+1)))
    for k in range(length):
        residual = (sum(f[i]*f[k-1-i] for i in range(k)) if k else 0)
        residual += sum(coeff(b, i)*f[k-i] for i in range(k+1))
        assert residual == (1 if k == 0 else 0)
    return f


def det(a):
    """Fraction-free Bareiss with row pivoting; includes singular and empty cases."""
    if not a:
        return 1
    a = [list(row) for row in a]
    size, old, sign = len(a), 1, 1
    for k in range(size-1):
        if not a[k][k]:
            pivot = next((i for i in range(k+1, size) if a[i][k]), None)
            if pivot is None:
                return 0
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        current = a[k][k]
        for i in range(k+1, size):
            for j in range(k+1, size):
                numerator = a[i][j]*current-a[i][k]*a[k][j]
                assert numerator % old == 0
                a[i][j] = numerator//old
            a[i][k] = 0
        old = current
    return sign*a[-1][-1]


def hankel(f, shift, size):
    return det([[f[shift+i+j] for j in range(size)] for i in range(size)])


def fidelity():
    examples = {
        2: [1, 1, -2, -1, 2, -1, -2, 1, 1, 0, 0, 0],
        3: [1, 0, -1, 2, 1, -1, 1, 0, -1, 0, 0, -1,
            0, 1, -1, 1, 2, -1, 0, 1, 0, 0, 0, 0],
    }
    count = 0
    for n, expected in examples.items():
        f = moments(n, n+2+4*len(expected))
        actual = [hankel(f, n+2, j) for j in range(2*len(expected))]
        assert actual[:len(expected)] == expected
        assert actual[len(expected):] == [(-1)**n*x for x in expected]
        count += len(actual)
        print('FIDELITY source Theorem B displayed block and next block:',
              'n=', n, 'block=', expected, flush=True)
    return count


def jet_checks():
    checks = 0
    for n in range(1, 9):
        data, u, r = blocks(n)
        rho = period(n)
        assert sum(k+1 for k, _, _ in data) == rho
        assert [[coeff(x, 0) for x in row] for row in r] == [[1, 0], [0, 1]]
        assert [[coeff(x, 0) for x in row] for row in u] == [[1, 0], [1, 0]]
        assert sub(mul(u[0][0], u[1][1]), mul(u[0][1], u[1][0])) == monomial(2*rho)
        cap = min(2*rho+3, 32)
        trace = add(u[0][0], u[1][1])
        # Hensel fixed point lambda_- = q^(2rho)/(trace-lambda_-).
        small = []
        for _ in range(cap):
            small = mul(monomial(2*rho), inv(sub(trace, small), cap), cap)
        large = sub(trim(trace[:cap]), small)
        gap_inv = inv(sub(large, small), cap)
        numerator = [[sub(u[i][j], small if i == j else []) for j in range(2)]
                     for i in range(2)]
        eplus = matscale(numerator, gap_inv, cap)
        eminus = [[sub(IDENTITY[i][j], eplus[i][j]) for j in range(2)] for i in range(2)]
        nil = [[sub(r[i][j], IDENTITY[i][j]) for j in range(2)] for i in range(2)]
        for k in (1, 2, 3, 7):
            rp = matpower(r, k, cap)
            binomial_sum = [[[], []], [[], []]]
            term = IDENTITY
            for a in range(cap):
                binomial_sum = matadd(binomial_sum,
                    matscale(term, [comb(k, a) if a <= k else 0], cap))
                term = matmul(term, nil, cap)
            assert truncate(rp, cap) == truncate(binomial_sum, cap)
            plus = matscale(eplus, power(large, k, cap), cap)
            minus = matscale(eminus, power(small, k, cap), cap)
            up = matpower(u, k, cap)
            assert truncate(up, cap) == truncate(matadd(plus, minus), cap)
            assert truncate(up, min(cap, 2*k*rho)) == truncate(plus, min(cap, 2*k*rho))
            checks += 3
        print('TRANSFER n=', n, 'rho=', rho, 'R(0)=I U(0)=[[1,0],[1,0]]',
              'det(U)=q^(2rho); finite jets OK', flush=True)
    return checks


def normal_inverse_checks():
    count = 0
    for n, size in [(1, 4), (1, 8), (2, 12), (2, 24)]:
        data = cycle(n)
        weight = sum(k+1 for k, _, _ in data)
        steps = len(data)*(size//weight)
        previous, current, degree = [], [1], 0
        for p in range(steps):
            k, v, d = data[p % len(data)]
            previous_k = data[(p-1) % len(data)][0]
            following = sub(mul(d, current),
                            mul(monomial(previous_k+k+2, v), previous))
            previous, current = current, following
            degree += k+1
        assert degree == size
        g = data[(steps-1) % len(data)][0]+1
        c = [coeff(current, size-i) for i in range(size+1)]
        d = [coeff(previous, size-g-i) for i in range(size)]
        f = moments(n, n+2+2*size+8)
        mu = f[n+1:]
        h = [[mu[i+j] for j in range(size)] for i in range(size)]
        normal_det = det(h)
        assert normal_det == (-1)**(n*(size//weight))
        assert all(sum(c[a]*mu[a+t] for a in range(size+1)) == 0 for t in range(size))
        assert all(sum(d[a]*mu[a+t] for a in range(size)) == (-1 if t == size-1 else 0)
                   for t in range(size))
        bezout = [[-sum(coeff(c, i+j+1-a)*coeff(d, a)-coeff(d, i+j+1-a)*coeff(c, a)
                       for a in range(j+1)) for j in range(size)] for i in range(size)]
        assert all(bezout[i][j] == bezout[j][i] for i in range(size) for j in range(size))
        assert all(sum(h[i][a]*bezout[a][j] for a in range(size)) == int(i == j)
                   for i in range(size) for j in range(size))
        zero_targets = 0
        for delta in range(2, 8):
            alpha, beta = delta//2, (delta+1)//2
            for j in range(size-beta+1):
                deleted_rows = [i for i in range(size) if not alpha <= i < alpha+j]
                deleted_cols = [i for i in range(size) if not beta <= i < beta+j]
                complement = det([[bezout[i][a] for a in deleted_rows] for i in deleted_cols])
                target = hankel(f, n+1+delta, j)
                assert target == (-1)**(j*delta)*normal_det*complement
                zero_targets += int(target == 0)
                count += 1
        print('NORMAL n=', n, 'N=', size, 'det=', normal_det, '; orthogonality, full Bezout inverse,',
              'all admissible complementary minors delta=2..7 OK;',
              'zero target minors=', zero_targets, flush=True)
    return count


def polynomial_checks():
    count = 0
    for n in (1, 2):
        rho = period(n)
        for delta, degree in ((2, 1), (3, 3)):
            shift = n+1+delta
            f = moments(n, shift+2*(6*rho-1)+1)
            max_degree = 0
            for residue in range(rho):
                values = [hankel(f, shift, residue+m*rho) for m in range(6)]
                differences, firsts = values[:], []
                while differences:
                    firsts.append(differences[0])
                    differences = [b-a for a, b in zip(differences, differences[1:])]
                assert all(x == 0 for x in firsts[degree+1:])
                # Newton data use m=0..degree; the remaining values are holdouts.
                assert all(values[m] == sum(firsts[a]*comb(m, a)
                           for a in range(degree+1)) for m in range(degree+1, 6))
                max_degree = max(max_degree, max((a for a, x in enumerate(firsts) if x), default=0))
                count += 6
            print('QUASIPOLYNOMIAL finite grid n=', n, 'ell=', shift, 'rho=', rho,
                  'all residues, m=0..5; Newton degree=', max_degree,
                  'holdouts per residue=', 5-degree, flush=True)
    return count


def main():
    print('Exact integer arithmetic; no floating-point evidence.', flush=True)
    print('Source: Han--Pedon arXiv:2610.36132v1, Theorem B and Conjecture C.', flush=True)
    source_count = fidelity()
    jet_count = jet_checks()
    complement_count = normal_inverse_checks()
    target_count = polynomial_checks()
    print('COUNTS source_determinants=', source_count, 'jet_identities=', jet_count,
          'complementary_minor_identities=', complement_count,
          'candidate_determinants=', target_count, flush=True)
    print('Finite corroboration only; universal Conjecture C remains unproved here.', flush=True)
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Exact finite checks for Hanna's OEIS A386667. No third-party dependencies.

Construct A from the source's A083563 reversion formula, independently of
any dyadic-support induction, and test the original functional equation
and literal quotient coefficients. All series are truncated modulo x^(M+1).
"""
from fractions import Fraction as Q

M = 128


def mul(a, b):
    out = [Q(0)] * (M + 1)
    aa = [(i, v) for i, v in enumerate(a) if v]
    bb = [(j, v) for j, v in enumerate(b) if v]
    for i, v in aa:
        for j, w in bb:
            if i + j > M:
                break
            out[i + j] += v * w
    return out


def power(a, n):
    out = [Q(1)] + [Q(0)] * M
    while n:
        if n % 2:
            out = mul(out, a)
        n //= 2
        if n:
            a = mul(a, a)
    return out


def inverse_unit(a):
    assert a[0] != 0
    b = [1 / a[0]] + [Q(0)] * M
    for n in range(1, M + 1):
        b[n] = -sum((a[j] * b[n-j] for j in range(1, n+1)), Q(0)) / a[0]
    return b


def compose(a, b):
    assert b[0] == 0
    out = [Q(0)] * (M + 1)
    p = [Q(1)] + [Q(0)] * M
    order = next((i for i, v in enumerate(b) if v), M + 1)
    for j in range(M // order + 1):
        if a[j]:
            out = [u + a[j] * v for u, v in zip(out, p)]
        p = mul(p, b)
    return out


def main():
    # Source (2.a): B=2x+(B^2+B(x^2))/2, A(x(1-B)^2)=x.
    B = [Q(0)] * (M + 1)
    for n in range(1, M + 1):
        B[n] = (4 * (n == 1) + (B[n//2] if n % 2 == 0 else 0)
                + sum((B[i]*B[n-i] for i in range(1, n)), Q(0))) / 2
    assert B[1:11] == list(map(Q, [2,3,6,18,54,183,636,2316,8610,32763]))
    V = [-b for b in B]
    V[0] = Q(1)
    U = power(inverse_unit(V), 2)
    A = [Q(0)] * (M + 1)
    p = [Q(1)] + [Q(0)] * M
    for n in range(1, M + 1):
        p = mul(p, U)
        A[n] = p[n-1] / n  # Exact formal reversion (Lagrange inversion).
    expected = [1,4,34,360,4271,54292,723030,9957288,140645083,2026330944]
    assert A[1:11] == list(map(Q, expected))
    print('SOURCE_FIDELITY A[1..10] =', ','.join(map(str, A[1:11])))
    print('SOURCE_FIDELITY A083563 B[1..10] =', ','.join(map(str, B[1:11])))
    f = [Q(0)] + power(V, 2)[:M]
    X = [Q(0), Q(1)] + [Q(0)] * (M-1)
    assert compose(A, f) == X
    assert compose(f, A) == X
    A2 = power(A, 2)
    inner = [x + 4*y for x, y in zip(X, A2)]
    assert compose(A, power(inner, 2)) == A2
    print('SOURCE_EQUATIONS original A^2=A((x+4A^2)^2), A(f)=x, f(A)=x: through degree', M)
    quotient_unit = A[1:] + [Q(0)]
    Phi = inverse_unit(quotient_unit)
    total = 0
    for n in range(2, 8):
        N = 2**n
        phiN = power(Phi, N)
        for k in range(n-1):
            a = 2**k
            value = phiN[N-a]
            inverse_value = power(f, a)[N]
            assert value * a == N * inverse_value
            assert value == 0, (n, k, value)
            assert inverse_value == 0, (n, k, inverse_value)
            print(f'DYADIC n={n} k={k} degree={N-a} coefficient={value} inverse_coefficient={inverse_value}')
            total += 1
    assert total == 21
    # Nonvacuity of the tested range: the immediately excluded first boundary is nonzero.
    boundary = power(Phi, 4)[2]
    assert boundary == 24
    print('BOUNDARY n=2 k=1 (excluded):', boundary)
    print('EXACT_CASES', total, '; no assertion beyond n=7 follows from this finite computation')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

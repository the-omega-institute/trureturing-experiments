"""Exact checks of the rank radius on cyclic difference sets and a collision."""

from itertools import combinations, product
from sympy import I, Matrix, Rational, simplify, sqrt, zeros


def min_excess(probabilities, rank):
    p = sorted(probabilities, reverse=True)
    assert sum(p) == 1 and min(p) >= 0 and rank >= 1
    if sum(1 for x in p if x > 0) <= rank:
        return 0
    for j in range(rank):
        level = sum(p[j:])/(rank-j)
        if p[j] <= level:
            return simplify((rank-j)*level**2-sum(x*x for x in p[j:]))
    raise AssertionError("No flattening index")


def verify_state(d, support, p, groups, expected):
    rho = zeros(len(p))
    for group in groups:
        v = Matrix([sqrt(w) if j in group else 0 for j, w in enumerate(p)])
        rho += v*v.H
    assert rho.diagonal() == Matrix([p]) and rho.rank() == len(groups)
    differences = [0]*d
    for j, k in product(range(len(p)), repeat=2):
        differences[(support[j]-support[k]) % d] += rho[j, k]
    squared = sum(simplify(v*v.conjugate()) for v in differences[1:])/d
    assert simplify(squared - expected/d) == 0
    assert min_excess(p, len(groups)) == expected


def concentrate(p, delta):
    result = sorted(p, reverse=True)
    remaining = min(delta, 1-result[0])
    result[0] += remaining
    for j in range(len(result)-1, 0, -1):
        amount = min(result[j], remaining)
        result[j] -= amount
        remaining -= amount
    return result


def check():
    for d, support in [(7, (0, 1, 3)), (13, (0, 1, 3, 9))]:
        assert len({(a-b) % d for a in support for b in support if a != b}) == len(support)*(len(support)-1)
    p = [Rational(1, 2), Rational(1, 3), Rational(1, 6)]
    verify_state(7, (0, 1, 3), p, [[0, 1, 2]], Rational(11, 18))
    verify_state(7, (0, 1, 3), p, [[0], [1, 2]], Rational(1, 9))
    p = [Rational(1, 2), Rational(1, 4), Rational(1, 8), Rational(1, 8)]
    for groups, excess in [([[0, 1, 2, 3]], 21), ([[0], [1, 2, 3]], 5), ([[0], [1], [2, 3]], 1)]:
        verify_state(13, (0, 1, 3, 9), p, groups, Rational(excess, 32))
    # Every diagonal of denominator eight in a fixed-support TV ball.
    target = [Rational(3, 8), Rational(2, 8), Rational(2, 8), Rational(1, 8)]
    lattice = []
    for i, j, k in combinations(range(11), 3):
        lattice.append([Rational(n, 8) for n in (i, j-i-1, k-j-1, 10-k)])
    for delta in (Rational(1, 8), Rational(2, 8)):
        ball = [p for p in lattice if sum(abs(a-b) for a, b in zip(p, target)) <= 2*delta]
        for rank in (1, 2, 3):
            assert min(min_excess(p, rank) for p in ball) == min_excess(concentrate(target, delta), rank)
    # At opposite sites the nonzero difference occurs twice, with opposite coherences.
    v = Matrix([sqrt(2)/2, 0, I*sqrt(2)/2, 0])
    for j in range(4):
        amplitude = sum(v[k]*I**(-j*k) for k in range(4))/2
        assert simplify(amplitude*amplitude.conjugate()) == Rational(1, 4)
    # In dimension three, leakage can cancel the inner support's two coherences.
    amplitudes = [Rational(2, 3), Rational(2, 3), Rational(-1, 3)]
    assert amplitudes[2]**2 == Rational(1, 9)
    assert amplitudes[0]**2 / (1-amplitudes[2]**2) == Rational(1, 2)
    for shift in range(3):
        correlation = sum(amplitudes[j] * amplitudes[(j-shift) % 3] for j in range(3))
        assert correlation == int(shift == 0)
    print("Rank witnesses, diagonal-noise thresholds, cyclic collision and leakage verified.")


if __name__ == "__main__":
    check()

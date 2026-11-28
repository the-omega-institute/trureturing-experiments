"""Exact finite differences and transport witnesses for noisy, tilted heat spectra."""

from fractions import Fraction
from math import comb, factorial, lcm


def moment(measure, order):
    return sum(weight * node**order for node, weight in measure.items())


def transport(left, right):
    nodes = sorted(left.keys() | right.keys())
    cumulative = 0
    cost = 0
    for a, b in zip(nodes, nodes[1:]):
        cumulative += left.get(a, 0) - right.get(a, 0)
        cost += abs(cumulative) * (b-a)
    return cost


def alternating(d, endpoint):
    m = 2*d-int(endpoint)
    step = Fraction(1, 4*m)
    start = 1-m*step if endpoint else Fraction(1, 4)
    mu, nu = {}, {}
    for j in range(m+1):
        (mu if j % 2 else nu)[start+j*step] = Fraction(comb(m, j), 2**(m-1))
    assert len(mu) == d and (1 in mu) == endpoint
    assert moment(mu, 0) == moment(nu, 0) == 1
    assert all(moment(mu, k) == moment(nu, k) for k in range(1, m))
    assert abs(moment(mu, m)-moment(nu, m)) == factorial(m)*step**m/2**(m-1)
    assert transport(mu, nu) == step
    # At u0=h, multiplicities are proportional to probability mass divided by x.
    weights = [mass/x for measure in (mu, nu) for x, mass in measure.items()]
    multiple = lcm(*(v.denominator for v in weights))
    assert all((v*multiple).denominator == 1 for v in weights)


def main():
    for d in range(1, 10):
        for endpoint in (False, True):
            alternating(d, endpoint)
    step = Fraction(1, 128)
    center = Fraction(1, 2)
    split = {center-step: Fraction(1, 2), center+step: Fraction(1, 2)}
    assert transport({center: 1}, split) == step
    assert moment(split, 1) == center and moment(split, 2)-center**2 == step**2
    for n in (2, 4, 8, 32):
        # In y units h=log(2), so x=1/2 is the mass coordinate y=1.
        eta = Fraction(1, n)
        tail = Fraction(1, 2**(n*n))
        mu = {Fraction(1, 2): 1}
        nu = {Fraction(1, 2): 1-eta, tail: eta}
        assert all(abs(moment(mu, j)-moment(nu, j)) <= eta for j in range(1, 19))
        assert transport({1: 1}, {1: 1-eta, n*n: eta}) == n-eta
    print("18 finite-difference spectra, separated-node splitting and unbounded mass transport verified.")


if __name__ == "__main__":
    main()

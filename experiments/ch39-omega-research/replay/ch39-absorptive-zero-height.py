"""Finite absorptive samples, confluent translated powers and threshold resolution."""
from math import isqrt, prod
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx, fmpq, fmpq_mat
from sympy import Rational, factor, integrate, symbols

FINITE = run_path(str(Path(__file__).with_name("ch39-finite-moment-height.py")))


def cut_kernels():
    beta, x = symbols("beta x")
    scalar = (1-2*x)**2/16
    vector = 3*scalar-Rational(1, 2)
    kernels = [scalar, vector, x*(1-x)]
    values = [integrate(kernel, (x, (1-beta)/2, (1+beta)/2)) for kernel in kernels]
    expected = [beta**3/48, beta**3/16-beta/2, beta/4-beta**3/12]
    assert all(factor(v-e) == 0 for v, e in zip(values, expected))


def falling(power, order):
    return prod(power-i for i in range(order))


def wronskians():
    # At s=100 all four square roots are rational. Compare every prefix of a
    # mixed simple/double basis to its closed confluent Vandermonde formula.
    s = 100
    nodes = [(19, 1), (36, 2), (51, 1), (64, 2)]
    columns = [(t, j) for t, mult in nodes for j in range(mult)]
    power = fmpq(3, 2)
    for k in range(1, len(columns)+1):
        pts = columns[:k]
        rows = []
        for i in range(k):
            row = []
            for t, derivative in pts:
                root = isqrt(s-t)
                assert root**2 == s-t
                entry = falling(power, i)*fmpq(root)**(3-2*i)
                if derivative:
                    entry *= -(power-i)/(s-t)
                row.append(entry)
            rows.append(row)
        expected = prod(falling(power, i) for i in range(k))
        for t, derivative in pts:
            expected *= isqrt(s-t)**3
            if derivative:
                expected *= fmpq(1, (s-t)**2)
        for j, (b, _) in enumerate(pts):
            for a, _ in pts[:j]:
                if a != b:
                    expected *= fmpq(1, s-b)-fmpq(1, s-a)
        assert fmpq_mat(rows).det() == expected and expected != 0


def resolution():
    ctx.prec = 180
    previous, limits = None, []
    epsilon = arb(342)/10**8
    for n in range(11485, 11491):
        gamma = acb.zeta_zero(n).imag
        ell = ((gamma/624)**2+arb(1)/4).sqrt()
        real = ((ell+arb(1)/2)/2).sqrt()
        imag = ((ell-arb(1)/2)/2).sqrt()
        D = acb(real, imag)
        E, W = D-1, -(arb(35)/13).sqrt()*(D+1)
        ratio = abs(E)**2/abs(W)**2
        y = (arb(35)/13*ratio).sqrt()
        inferred = 2*y/(1-y)**2
        height = 624*(inferred**2-arb(1)/4).sqrt()
        assert height.overlaps(gamma)
        if previous is not None:
            u = (ratio/previous).sqrt()
            limits.append((u-1)/(u+1))
            assert previous*(1+epsilon)/(1-epsilon) < ratio*(1-epsilon)/(1+epsilon)
        print(n, "extreme fermion threshold ratio", ratio.str(15))
        previous = ratio
    assert all(limits[-1] < v for v in limits[:-1])
    assert arb("3.42545090") < limits[-1]*10**6 < arb("3.42545092")
    print("Worst independent relative threshold error, boundary in ppm",
          (limits[-1]*10**6).str(13))


if __name__ == "__main__":
    cut_kernels()
    wronskians()
    factors = FINITE["charged_scalar_factors"]()
    count = sum(FINITE["Poly"](p, FINITE["x"]).degree() for p in factors)+3
    assert 2*(count+2*(3+4)) == 52
    resolution()

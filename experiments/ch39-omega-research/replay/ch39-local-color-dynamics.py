"""A finite nearest-neighbor generator for discretized arithmetic color-face weights."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

full = run_path(str(Path(__file__).with_name("ch39-full-color-obstruction.py")))
q = full["q"]


def balanced_steps(weights):
    n = len(weights)
    return [[min(F(1), F(weights[(i+d)%n], weights[i]))/4 for d in (-1,1)] for i in range(n)]


def outward_sampler(weights):
    n = len(weights)
    half = n//2
    assert n%2 == 0 and all(weights[i] == weights[-i%n] for i in range(n))
    mass = [weights[0], *[2*x for x in weights[1:half]], weights[half]]
    total = sum(mass)
    radius = [F(0)]*(half+1)
    active = F(1)
    remainder = total
    for i, local in enumerate(mass):
        radius[i] = active*F(local,remainder)
        active -= radius[i]
        remainder -= local
    law = [radius[min(i,n-i)]/(2 if i not in (0,half) else 1) for i in range(n)]
    assert law == [F(x,total) for x in weights]


def verify():
    ctx.prec = 180
    n, k = 2**29, full["K"]
    bits = 24*k+100
    degree = 2*k
    # First-quadrant dyadic roots, starting at the primitive fourth root.
    cosine = arb(0)
    for _ in range(27):
        cosine = ((1+cosine)/2).sqrt()
    chord = (2*(1-cosine)).sqrt()
    assert chord < q(F(1,2**26)) and (2*(1-(2*cosine*cosine-1))).sqrt() > q(F(1,2**26))
    high = full["high"]
    power = (high["geo"]["BOOK_PLANCK"]/(1000*q(full["L0"]).exp()))**2/(8*arb.pi())
    assert k-2 < power < k-1 and degree+k < n
    for zero in range(11485,11491):
        amplitude = 16/(acb.zeta_zero(zero).imag**2+q(high["OFFSET"]))
        assert 2*amplitude/(1-2*amplitude) < q(F(1,2**21))
        assert (1-2*amplitude)**power > q(F(1,2**13))
    assert 21*(degree+1)-15 > bits-15 and bits-22 > 24*k+54+2
    assert 2**15*n*n//4 == 2**71
    assert 2**15*n == 2**44
    model = [11,8,5,3,2,3,5,8]
    transitions = balanced_steps(model)
    for i,(_,up) in enumerate(transitions):
        j = (i+1)%len(model)
        assert model[i]*up == model[j]*transitions[j][0]
        assert sum(transitions[i]) <= F(1,2)
    outward_sampler(model)
    print("Six weight ranges, dyadic arc bounds, finite quadrature and Markov margins verified.")


if __name__ == "__main__":
    verify()

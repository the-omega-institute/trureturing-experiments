"""A finite unitary feedback channel for the threshold phase of the same shell zero."""
from fractions import Fraction as F
from math import gcd,isqrt
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

ROOT = Path(__file__).parent
PHASE = run_path(str(ROOT/"ch39-balanced-mass-phase.py"))
a = PHASE["a"]


def certify_masses(phase, r, z):
    d,h = [(a(x)*z).sqrt()*r for x in ("13/10","7/2")]
    t = d/abs(phase-1).sqrt()
    real = ((t*t+d*d/2)/2).sqrt()
    imag = (t*t-real*real).sqrt()
    for sign in (1,-1):
        middle = sign*acb(real,imag)
        c,e = middle+d,middle-d
        p,w = 1/(abs(c)*t*t),1/t**3
        if not (abs(p-w)<h<p+w):
            continue
        wr = (p*p-w*w-h*h)/(2*h)
        weak = acb(wr,(w*w-wr*wr).sqrt())
        color = weak+h
        assert (c*e/middle**2-phase).contains(0)
        assert (abs(color*c*middle**2)-1).contains(0)
        assert (abs(weak*middle**3)-1).contains(0)
        return t
    raise AssertionError("No mass intersection")


def verify():
    ctx.prec = 150
    # Primitive rational symmetric orthogonal couplers up through denominator 85.
    candidates = []
    for n in range(1,86):
        for u in range(1,n):
            v = isqrt(n*n-u*u)
            if u*u+v*v==n*n and gcd(u,v)==1 and F(n+u,n-u)>131:
                candidates.append((u,v,n))
    assert candidates==[(84,13,85)]
    smallest_rational = F(84,85)
    shell_hop = F(1,19)
    shell_coupling = (1-shell_hop**2)/(1+shell_hop**2)
    assert shell_coupling == F(180,181)
    for n,kind,r,v,state in PHASE["MATCH"]["matched_points"]():
        zero = acb.zeta_zero(n)
        gamma = zero.imag
        w = (zero-1)/zero
        if n==11490:
            _,min_chord,_ = PHASE["phase_bounds"](r,1/state[0])
            assert 1/((gamma/131)**2+a("1/4")).sqrt() < min_chord
        for refl in (F(131,133),smallest_rational,shell_coupling):
            gain = (1+refl)/(1-refl)
            reflection = a(refl)
            transmission = a(1-refl*refl).sqrt()
            inner = transmission*w/(1-reflection*w)
            output = transmission*inner-reflection
            reduced = acb(a("1/2"),gamma/a(gain))
            assert (output-(reduced-1)/reduced).contains(0)
            assert (abs(output)-1).contains(0)
            t = certify_masses(output,r,1/state[0])
            if refl==shell_coupling:
                assert gain==361 and (transmission-a("19/181")).contains(0)
                assert a("18.9") < abs(inner) < 19
                print(n,kind,"effective height",(gamma/361).str(11),
                      "output chord",abs(output-1).str(9),"|D|",t.str(8))


if __name__=="__main__":
    verify()

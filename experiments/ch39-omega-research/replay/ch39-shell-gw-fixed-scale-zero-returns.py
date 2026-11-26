"""All near returns at fixed one-loop mass on the critical line up to a finite height."""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

fixed = run_path(str(Path(__file__).with_name("ch39-shell-gw-fixed-scale-density.py")))
q = arb


def candidates(f, ratio, center, cutoff, width, radial):
    period = 2*q.pi()
    # (denominator, numerator) is a lattice basis; its determinant is exactly ±1.
    q1, p1 = 335627933, 614345469
    q2, p2 = 242126276, 443196665
    assert abs(p1*q2-p2*q1) == 1
    e1, e2 = q1*ratio-p1, q2*ratio-p2
    limit = int(((f*(cutoff-center)+width)/period).ceil().unique_fmpz())
    tol = (ratio*width+radial)/period
    def ceiling(v):
        return int(v.upper().ceil().unique_fmpz())
    ca, cb = ceiling(limit*abs(e2)+q2*tol), ceiling(limit*abs(e1)+q1*tol)
    hits = []
    for a in range(-ca, ca+1):
        for b in range(-cb, cb+1):
            k, m = a*q1+b*q2, a*p1+b*p2
            if 0 < k <= limit and not abs(a*e1+b*e2)>tol:
                hits.append((k,m,a*e1+b*e2))
    return sorted(hits)


def exclude(s, half, depth=0):
    z = acb(q(1)/2,s+arb(0,half.upper())).zeta()
    if abs(z).abs_lower()>0:
        return abs(z).abs_lower()
    assert depth<8
    return min(exclude(s-half/2,half/2,depth+1),exclude(s+half/2,half/2,depth+1))


def certify(index, cutoff):
    ctx.prec = 128
    a = [(1-q(d)/26).acos() for d in (3,5,7)]
    f, g = [angle-a[0] for angle in a[1:]]
    ratio = g/f
    center = acb.zeta_zero(index).imag
    assert acb.zeta_zero(index+1).imag-center>q(".01")
    assert center-acb.zeta_zero(index-1).imag>q(".01")
    if index==43:
        width, lower, tilt = q("3e-8"), q("32.63"), q(".007")
    else:
        width, lower, tilt = q("1e-8"), q("12.816"), q(".0035")
    exact = tilt*width/lower
    hits = candidates(f,ratio,center,cutoff,width,exact)
    k,m,error = hits[0]
    assert 2*q.pi()*abs(error)<ratio*width-exact
    first = center+2*q.pi()*m/g + arb(0,(exact/g).upper())
    print("window",index,"first continuous return",k,m,first.str(24))
    # Count the narrower tolerance strips, then exclude zeros on the full square.
    radial = (tilt*width+q("1e-8"))/lower
    narrow = candidates(f,ratio,center,cutoff,width,radial)
    hits = candidates(f,ratio,center,cutoff,width,width)
    assert center-(2*q.pi()*hits[0][0]-width)/f < 0
    minimum = q(100)
    for _,m,_ in hits:
        s = center+2*q.pi()*m/g
        minimum = min(minimum,exclude(s,width/g))
    print("exposure tolerance 1e-8 and full-window candidate counts",len(narrow),len(hits),
          "whole-window |zeta| lower on off-center strips",minimum.str(7))
    assert minimum > q(".009" if index == 43 else ".155")


if __name__ == "__main__":
    for n in (43,195):
        certify(n,q("1e11"))

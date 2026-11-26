"""Reflected triangle returns, a reciprocal one-port, and two-port probabilities."""
from collections import Counter
from itertools import product

from flint import acb, arb, ctx


def check_orbits():
    pts = {p for p in product(range(30),repeat=3) if sum(x*x for x in p)==867}
    def rot(p):
        return p[1],p[2],p[0]
    near = (13,13,23)
    marked = [near,rot(near),rot(rot(near))]
    ring = [(marked[n%3],n//3) for n in range(9)]
    tags = {(p,j) for p in pts for j in range(3 if p in marked else 1)}
    def g(u):
        return ring[(ring.index(u)+1)%9] if u in ring else (rot(u[0]),0)
    def refl(u):
        p,_ = u
        return ring[-ring.index(u)%9] if u in ring else ((p[1],p[0],p[2]),0)
    assert all(refl(refl(u))==u and g(refl(g(u)))==refl(u) for u in tags)
    p = ((17,17,17),0)
    todo = set(product(tags,repeat=2))-{(p,p)}
    stable = Counter()
    while todo:
        u = min(todo)
        cycle = []
        while u in todo:
            todo.remove(u)
            cycle.append(u)
            u = tuple(map(g,u))
        if tuple(map(refl,u)) in cycle:
            stable[len(cycle)] += 1
    assert stable=={3:3,9:5}


def check_probabilities():
    ctx.prec = 150
    probs = []
    for n in range(11485,11491):
        zero = acb.zeta_zero(n)
        gamma = zero.imag
        w = (zero-1)/zero
        gain = 624
        kappa = arb(gain-1)/(gain+1)
        b = (w-kappa)/(1-kappa*w)
        conversion = 16*gamma**2*gain**2/(4*gamma**2+gain**2)**2
        assert (b.imag**2-conversion).contains(0)
        assert (b.real**2+b.imag**2-1).contains(0)
        probs.append(conversion)
    gaps = [x-y for x,y in zip(probs,probs[1:])]
    assert all(val>arb("2.5995e-7") for val in gaps) and gaps[-1]<arb("2.5996e-7")
    assert arb(".00313342288") < probs[-1] < arb(".00313342289")
    print("Conversion probabilities", [p.str(12) for p in probs],
          "minimum gap",gaps[-1].str(10))
    assert all(gaps[-1]<g for g in gaps[:-1])


if __name__ == "__main__":
    check_orbits()
    check_probabilities()

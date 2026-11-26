"""Local reflection, clocked feedback and a finite time-bin mirror readout."""
from collections import deque
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx
from sympy import Rational, cancel, symbols

LOCAL = run_path(str(Path(__file__).with_name("ch39-local-shell-return.py")))


def distances(graph, start):
    todo, dist = deque([start]), {start:0}
    while todo:
        p = todo.popleft()
        for q in graph[p]:
            if q not in dist:
                dist[q] = dist[p]+1
                todo.append(q)
    return dist


def radii(den, cutoff):
    pole = (17*den,)*3
    branches = LOCAL["graph"](den,cutoff)
    endpoints = {LOCAL["scaled"](p,den) for p in LOCAL["shell_grid"](1)} | {pole}
    sites = set(branches) | {pole}
    for flat in (False,True):
        graph = {p:[q for q,w in branches.get(p,[]) if not flat or w==0] for p in sites}
        for p in sites-{pole}:
            if LOCAL["squared"](p,pole)<=cutoff*den**2:
                graph[p].append(pole)
                graph[pole].append(p)
        ecc = {}
        for p in sites:
            distance = distances(graph,p)
            ecc[p] = max(distance[q] for q in endpoints)
        radius = min(ecc.values())
        assert [p for p in sites if ecc[p]==radius] == [pole]
        assert (radius,max(ecc[q] for q in endpoints)) == ((4,6) if den==3 else (2,4))
    return radius


def verify():
    assert radii(3,56)==4 and radii(1,152)==2
    z, w = symbols("z w")
    d = Rational(623,625)
    ret = z**17*w
    amplitude = z**8*(-d+(1-d**2)*ret/(1-d*ret))
    assert cancel(amplitude-z**8*(ret-d)/(1-d*ret)) == 0
    ctx.prec = 180
    d = arb(623)/625
    m = 67**2-1
    assert 8+17*m==76304
    steady, finite = [], []
    for n in range(11485,11491):
        rho = acb.zeta_zero(n)
        phase = (rho-1)/rho
        direct = (phase-d)/(1-d*phase)
        b = direct-(1-d**2)*phase*(d*phase)**m/(1-d*phase)
        # Independent finite input recurrence, restricted to the uniform internal mode.
        echo_sum = acb(0)
        for _ in range(m):
            echo_sum = 1+d*phase*echo_sum
        assert b.overlaps(-d+(1-d**2)*phase*echo_sum)
        assert abs(b-direct) < 2*d**m
        steady.append(direct.imag**2)
        finite.append(b.imag**2)
        assert abs(finite[-1]-steady[-1]) < arb("1.28e-7")
        delay = 8+17*(1-d**2)/abs(1-d*phase)**2
        assert arb("10607.67") < delay < arb("10607.69")
    assert all(steady[i]-steady[i+1] > arb("2.59e-7") for i in range(5))
    gaps = [u-v for u,v in zip(finite,finite[1:])]
    assert all(g > arb("2.59e-7") for g in gaps)
    print("Designated output bin",8+17*m,"incident bins",m+1)
    print("Finite calibrated intensity gaps",[g.str(10) for g in gaps])
    print("Smallest unconditional time-bin probability gap",
          min(g.lower() for g in gaps)/(m+1))


if __name__ == "__main__":
    verify()

"""Derive arithmetic doublet weights from the maximal-return fibers of the 624-channel shell."""
from fractions import Fraction
from itertools import product
from pathlib import Path
from runpy import run_path
from flint import acb, arb, ctx

base = Path(__file__).parent


def fiber_weights():
    shell = run_path(str(base/"ch39-rotor-shell-returns.py"))
    nodes, adj, _ = shell["charged_shell"]()
    amplitudes = [shell["returns"](adj,i)[2].get(1,Fraction()) for i in range(len(nodes))]
    selected = {tuple(x//3 for x in p) for p,b in zip(nodes,amplitudes) if b==max(amplitudes)}
    assert len(selected)==3
    geometry = run_path(str(base/"ch39-intrinsic-shell-comparisons.py"))
    distance = geometry["squared"]
    points = {p for p in product(range(30),repeat=3) if sum(x*x for x in p)==867}
    pole = (17,17,17)
    nearest = {v for v in points if distance(v,pole)==68}
    assert len(nearest)==3 and not selected & nearest
    total, marked = 0, []
    for p in points-{pole}:
        values = {sum((x-17)*(y-17) for x,y in zip(p,r)) for r in nearest}
        signed = values | {-v for v in values}
        size = len(signed)**2
        total += size*(3 if p in nearest else 1)
        if p in selected:
            assert signed == {-96,-44,44,96}
            marked.append(size)
    assert total==624 and marked==[16]*3
    return Fraction(len(marked),total), Fraction(sum(marked),total)


def verify():
    ctx.prec = 240
    trace, coherent = fiber_weights()
    assert (trace,coherent)==(Fraction(1,208),Fraction(1,13))
    assert 2/Fraction(624)<Fraction(".00364249896012")<Fraction(".00364249896016")<trace
    finite = run_path(str(base/"ch39-finite-arithmetic-height.py"))
    g = acb.zeta_zero(10875).imag
    best = finite["weight"](g)
    maximum = finite["source"]["modulus"](-arb(59)/108)/(arb(10603)**2+arb(1)/4)
    h = arb(coherent.numerator)/coherent.denominator
    c = arb(trace.numerator)/trace.denominator
    assert maximum/best < arb("1.003451")
    assert maximum/best*(1-h/2) < 1-arb("1e-10")
    assert h/(2-h)>c>arb(1)/256
    grid = 2**64
    numerator = int((grid*g).floor().unique_fmpz())
    rows = []
    for m in range(numerator-10,numerator+12):
        y = arb(m)/grid
        t = abs(finite["finite_value"](y))
        w = finite["weight"](y)
        rows.append((m,w/(1+c*t),w*(1-h+h/(1+t))))
    for col, coefficient in [(1,c),(2,h/(2-h))]:
        winner = max(rows,key=lambda row:row[col].mid())
        assert winner[0]==195923358475422178854846
        assert all(winner[col]>row[col] for row in rows if row[0]!=winner[0])
        l = winner[col].lower()/best
        assert l>1-arb("1e-10")
        radius = (1-l+l*coefficient*arb("3.65e-22"))/(5*l*coefficient-arb(".015"))
        assert 0 < radius < 10/arb(grid)
    print("Maximal-return fibers",3,"normalized trace",trace,"pure equal-channel weight",coherent)
    print("Both finite responses have unique 64-bit height numerator",winner[0])


if __name__ == "__main__":
    verify()

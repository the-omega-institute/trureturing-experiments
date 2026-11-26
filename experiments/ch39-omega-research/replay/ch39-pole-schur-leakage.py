"""A nonreal three-clock Schur pivot detects additional mass paths at the external pole."""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx

REFINE = run_path(str(Path(__file__).with_name("ch39-refined-mass-obstruction.py")))
ball = REFINE["ball"]


def verify():
    ctx.prec = 170
    nodes, edges = REFINE["geometry"]()
    pole = nodes.index((51,51,51))
    incident = [(j if i==pole else i, weight) for i,j,weight,_ in edges if pole in (i,j)]
    signature = Counter(REFINE["LOCAL"]["squared"](nodes[pole],nodes[j]) for j,_ in incident)
    assert signature == {204:6,408:6}
    total = sum(weight**2 for _,weight in incident)
    assert total == F(103630995,278837197088)
    r = ball(F(98,1000))
    eta = arb(3).sqrt()/2
    clock = acb(-arb(1)/2, eta)
    lower = eta*ball(total)/((r+arb(1)/2)**2+eta**2)
    upper = ball(total)/eta
    assert lower > arb(".00029059") and upper < arb(".00043")
    assert ball(total)**2/(1+r**3) > arb("1.3799e-7")
    number = len(nodes)
    eye = acb_mat([[int(i==j) for j in range(number)] for i in range(number)])
    for n in range(11485,11491):
        rho = acb.zeta_zero(n)
        M = REFINE["matrices"](nodes,edges,(rho-1)/rho)
        full = (clock*eye-M).inv()
        column_norm = sum(abs(full[j,pole])**2 for j in range(number) if j!=pole)
        excess = (1/full[pole,pole]).imag-eta
        assert lower < excess < upper
        assert (excess-eta*column_norm/abs(full[pole,pole])**2).contains(0)
        mass = M*(eye-M**3).inv()
        assert mass[pole,pole].real > arb("1.3799e-7")
        print(n,"imaginary pivot excess",excess.str(14),"pole mass response",mass[pole,pole].real.str(14))


if __name__ == "__main__":
    verify()

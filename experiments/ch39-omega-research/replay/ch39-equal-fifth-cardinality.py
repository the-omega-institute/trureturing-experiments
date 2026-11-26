"""The exact fifty-six-state minimum and a shell-symmetric realization on three sheets."""
from fractions import Fraction as F
from itertools import product
from collections import Counter
from pathlib import Path
from runpy import run_path
from flint import arb, acb, arb_mat, ctx

ROOT = Path(__file__).resolve().parent
fifth = run_path(str(ROOT/"ch39-fifth-shell-harmonic.py"))
two, q = fifth["two"], fifth["q"]
SEEDS = {
    11485: (5268184311758, 8266506928826, 8727121110202, 9004139020976),
    11486: (5267098881711, 8266193794432, 8728849996082, 9003998090659),
    11487: (5265441304960, 8265716041586, 8731493253194, 9003782503954),
    11488: (5264727877133, 8265510577050, 8732632155024, 9003689574557),
    11489: (5263220512670, 8265076780124, 8735041144266, 9003492946350),
    11490: (5262508224967, 8264871943813, 8736180822714, 9003399896509),
}


def extremal_bound(center, moments):
    weights = [q(F(k, 55)) for k in (1,24,1,29)]
    powers = (1,4,9,16)
    radius = q("1e-8")
    box = [x+arb(0,radius) for x in center]
    assert q(0) < box[0] < box[1] < box[2] < box[3] < 1
    def derivative(nodes):
        return arb_mat([[w*k*v**(k-1) for w,v in zip(weights,nodes)] for k in powers])
    middle = derivative(center).inv()
    pre = arb_mat(4,4,[v.mid() for v in middle.entries()])
    assert not pre.det().contains(0)
    jac = arb_mat([[int(i==j) for j in range(4)] for i in range(4)])-pre*derivative(box)
    residual = arb_mat([[sum(w*x**k for w,x in zip(weights,center))-target]
                        for k,target in zip(powers,moments[1:5])])
    delta = pre*residual
    for i in range(4):
        assert abs(delta[i,0])+radius*sum(abs(jac[i,j]) for j in range(4)) < radius
    return sum(w*x**25 for w,x in zip(weights,box))-moments[5]


def shell_assignment():
    shell = [p for p in product(range(30), repeat=3) if sum(x*x for x in p)==867]
    # For each sheet read pole, inner orbit, and complete fifteen-point outer shell.
    sheets = [(1,3,3),(2,3,5),(4,4,5)]
    row = {51:0,49:1,47:2,35:2}
    result = Counter(table[row[sum(p)]] for p in shell for table in sheets)
    assert [result[i] for i in range(1,6)] == [1,1,21,4,30]
    assert sum(y==z and x+y+z in (35,47) for x,y,z in shell)==1


def sign_changes():
    # A repeated reference value allows a negative block followed by a positive one.
    single = [(),(-1,),(1,)]
    repeated = [*single,(-1,1)]
    for pieces in product(single, repeated, single, repeated):
        signs = sum(pieces, ())
        changes = sum(a!=b for a,b in zip(signs,signs[1:]))
        assert changes <= 4
        if changes == 4:
            assert signs[0] == signs[-1] == 1


def verify():
    ctx.prec = 250
    sign_changes()
    shell_assignment()
    mu = 1000*q(F(429,16)).exp()
    power = (two["geo"]["BOOK_PLANCK"]/mu)**2/(8*arb.pi())
    for index, seed in SEEDS.items():
        goal = two["harmonics"](acb.zeta_zero(index).imag,power,155,5)
        boxes = fifth["moment_boxes"](goal,fifth["POWERS"])
        margin = extremal_bound([q(F(x,10**13)) for x in seed],boxes)
        assert margin > q("0.00000160")
        print("Zero, fifth ideal moment separation at 55 equal labels:", index, margin.str(9))


if __name__ == "__main__":
    verify()

"""Exclude uniform site potentials for a rectangle of radial kinetic weights."""
from itertools import product
from pathlib import Path
from runpy import run_path
from flint import acb, acb_mat, arb

folder = Path(__file__).resolve().parent
common = run_path(str(folder/"ch39-shell-gw-195-common-field.py"))
first = run_path(str(folder/"ch39-shell-gw-native-kinetic-cost.py"))
model = common["flavor"]
adj, rows, count = common["adj"], common["rows"], common["count"]


def radial_form(kinetic, weights):
    points = model["pts"]
    labels = [[sum(a*b for a,b in zip(points[i],points[j])) for j in range(6)]
              for i in range(6)]
    degrees = [sum(weights.get(k,0) for k in row) for row in labels]
    form = acb_mat(54,54)
    for i,j in product(range(54),repeat=2):
        if i == j:
            form[i,j] = degrees[i//9]
        elif i//9 != j//9:
            form[i,j] = weights.get(labels[i//9][j//9],0)*kinetic[i,j]
    return form


def gap(index):
    if index==43:
        K, F = first["scalar_kinetic"](), first["m"]["F"]
    else:
        K, F = common["operators"]()
    t = common["native"]["relative_phase"](arb(676))
    phase = [acb(0,(-1)**i*t/2).exp() for i in range(6)]
    copy = acb_mat([[phase[i//9] if i%9==j else 0 for j in range(9)] for i in range(54)])
    G = adj(F)*F
    Y = common["conjugate"](F,t)
    D = adj(Y)*Y
    u2,d2 = common["old"]["u2"],common["old"]["d2"]
    x,b = arb(350),arb(20000)
    pencil = b*(D-d2/u2*G)-x*G
    def test(weights,bound,pencil0,level):
        form = radial_form(K,weights)
        block = adj(copy)*form*copy/6
        for i in range(54):
            for j in range(9):
                form[i,9*(i//9)+j] -= block[i%9,j]
        assert count(rows(form+pencil0),bound-level*u2/(4*arb.pi())**2)==0
    bound = arb(".100") if index==43 else arb(".067")
    for w19,w22 in product((arb(".13"),arb(1)),(arb(".24"),arb(1))):
        test({19:w19,22:w22,25:1},bound,pencil,x)
    x = arb(500)
    bound_equal = arb(".972") if index==43 else arb(".966")
    test({19:1,22:1,25:1},bound_equal,30000*(D-d2/u2*G)-x*G,x)
    print("zero",index,"uniform excess above",bound,"equal-edge excess above",bound_equal)


if __name__ == "__main__":
    theta=lambda h:(arb(h)/26).acos()
    for h,floor in ((19,".13"),(22,".24")):
        assert ((theta(25)/theta(h))**2)>arb(floor)
        assert arb(1)/(26-h)>arb(floor)
    for index in (43,195):
        gap(index)

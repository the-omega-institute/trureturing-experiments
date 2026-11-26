"""Compute oriented flavor gain using the original shell's two clocked portals."""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb

root = Path(__file__).parent
fit = run_path(str(root/"ch39-static-flavor-action.py"))
geo = run_path(str(root/"ch39-arc-relative-determinant.py"))
axes = fit["shell"]["axes"]
c = fit["c2"].sqrt()
origin, first, second = (geo["sites"].index(p) for p in
                         ((1,5,29), (13,23,13), (23,13,13)))
inside = [i for i in range(len(geo["sites"])) if i != origin]
source = acb_mat([[int(i == second)] for i in inside])


def signed_branches():
    n, m = axes
    norm_product = arb(sum(x*x for x in n)*sum(y*y for y in m)).sqrt()
    comm13 = c*(n[0]*m[2]-n[2]*m[0])/norm_product
    near, next_opposite = [], []
    for xy in fit["enclose_roots"]():
        w = fit["rephased_matrix"](xy)
        tr = sum(w[i][i] for i in range(3)) - 1
        hd = (abs(tr)/2).acos()
        skew = ((w[0][2]-w[2][0])*tr.conjugate()/abs(tr)).real
        coefficient = -hd*skew/(2*hd.sin()*comm13)
        near.append(coefficient)
        next_opposite.append(-coefficient*(arb.pi()/hd-1))
    assert all(abs(l) > arb(".12") for l in next_opposite)
    assert near[3] < 0 and near[0] < near[3] and 0 < near[2] and near[2] < near[1]
    return near


def shell_coefficient(height):
    matrix = geo["matrix"](0, geo["phases"](height))
    block = acb_mat([[matrix[i,j] for j in inside] for i in inside])
    response = block.solve(source)
    exchange = -matrix[origin,first]*response[inside.index(first),0]*matrix[second,origin]
    return exchange.imag/c


def min_gain(value, branches):
    if value > 0:
        return branches[2]/value
    assert value < 0
    return branches[3]/value


def verify():
    near = signed_branches()
    gains = {}
    for index, zero in enumerate(acb.zeta_zeros(10815,1278),10815):
        value = shell_coefficient(zero.imag)
        assert abs(value) < arb(".054")
        assert all(0 not in value-l for l in near)
        gains[index] = min_gain(value, near)
        if index in (11439,11737):
            assert abs(value) < arb(".03010")
            print(index, "coefficient", value.str(12), "least gain", gains[index].str(12))
    assert arb("3.18344") < gains[11439] and gains[11439] < arb("3.18345")
    assert arb("10.58355") < gains[11737] and gains[11737] < arb("10.58356")
    winner = gains[11431]
    assert arb(".72539") < winner and winner < arb(".72540")
    assert all(winner < value for idx,value in gains.items() if idx != 11431)
    print("unique least gain among 1278 candidates: 11431", winner.str(12))


if __name__ == "__main__":
    verify()

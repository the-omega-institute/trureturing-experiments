"""Directional color readout from scalar energy and a robust 56-state lower bound."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx, fmpq, fmpq_mat

root = Path(__file__).with_name
scalar = run_path(str(root("ch39-scalar-block-metric.py")))
equal = run_path(str(root("ch39-equal-fifth-cardinality.py")))
fifth,q,two = equal["fifth"],equal["q"],equal["two"]


def scalar_checks():
    a,b = fmpq(283,540),-fmpq(11,54)
    Q = fmpq_mat([[1+int(i==j) for j in range(4)] for i in range(4)])
    e = fmpq_mat([[1,0,0,0]])
    total = fmpq_mat(4,4)
    for order in scalar["orders"]:
        det,_,inv = scalar["simplex"](order,0)
        assert det == fmpq(1,4)
        total += inv/48
    assert total == a*Q+b*e.transpose()*e
    assert 3*b/25 == -fmpq(11,450)
    t = (0,1,4,9,36)
    table = fmpq_mat([[x**j for x in t] for j in range(5)])
    assert table.det() != 0


def robust_cardinality():
    ctx.prec = 320
    scale = 1000*q(F(429,16)).exp()
    power = (two["geo"]["BOOK_PLANCK"]/scale)**2/(8*arb.pi())
    nodes = [(q(".827"),False),(q(".827"),True),(q(".9"),False),(q(".9"),True)]
    cut = q(".527")
    low = fifth["constrained_polynomial"]([(cut,False),*nodes],fifth["POWERS"],-1)
    cap = fifth["constrained_polynomial"](nodes,fifth["POWERS"][:-1],1)
    at_cut = sum(a*cut**k for a,k in zip(cap,fifth["POWERS"]))
    for index,seed in equal["SEEDS"].items():
        target = two["harmonics"](acb.zeta_zero(index).imag,power,170,5)
        boxes = fifth["moment_boxes"](target,fifth["POWERS"])
        boxes = [boxes[0],*[m+arb(0,q("1e-12")) for m in boxes[1:]]]
        assert sum(c*m for c,m in zip(low,boxes)) > 0
        assert sum(c*m for c,m in zip(cap,boxes))/at_cut < q(F(1,54))
        center = [q(F(v,10**13)) for v in seed]
        assert equal["extremal_bound"](center,boxes) > q(".00000160")
    print("Six five-harmonic boxes of radius 1e-12 exclude every uniform size below 56.")


if __name__ == "__main__":
    scalar_checks()
    robust_cardinality()

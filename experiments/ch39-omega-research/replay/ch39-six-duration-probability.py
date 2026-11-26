"""保留全部原味系数，用正时长和二值跃迁概率识别纯 11431。"""
from fractions import Fraction
from math import factorial, prod
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx

certificate = run_path(str(Path(__file__).with_name("ch39-exposed-flavor-response.py")))
ordered = run_path(str(Path(__file__).with_name("ch39-ordered-flavor-unitary.py")))
study = certificate["study"]
labels = [11381,11397,11430,11432,11303,11576]
ticks = [2524400588,27079266863,228067471734,631469107462,54654922957,56204830396]
assert sum(ticks) == 10**12


def ball(q):
    return arb(q.numerator)/q.denominator


def cross(x, y):
    return [x[(i+1)%3]*y[(i+2)%3]-x[(i+2)%3]*y[(i+1)%3] for i in range(3)]


def projector(axis):
    x, y, z = axis
    return acb_mat([[1+z, acb(x,-y)], [acb(x,y), 1-z]])/2


def verify():
    ctx.prec = 192
    A, B, _ = certificate["polynomial"]()
    normal = [B, arb(1), A]
    length = sum(v*v for v in normal).sqrt()
    r = [1/(1+B*B).sqrt(), -B/(1+B*B).sqrt(), arb(0)]
    m = cross(r, [v/length for v in normal])
    rho, event = projector(r), projector(m)
    for unit in ([1,0,0], [0,1,0], [0,0,1]):
        pauli = 2*projector(unit)-ordered["flavor"]["identity"]
        slope = (acb(0,1)*event*(pauli*rho-rho*pauli)).trace().real
        expected = sum(u*v for u,v in zip(unit, normal))/length
        assert abs(slope-expected) < arb("1e-32")

    eta, h, k, _ = study["coefficients"]
    P, Q = [ordered["flavor"]["arcs"][d] for d in (32,152)]
    overlap = (P*Q).trace().real
    shell_bound = arb(9)/4
    K = arb(11)/4
    assert eta > 0 and (eta*shell_bound +
            (h*h+k*k+2*h*k*(2*overlap-1)).sqrt()) < K
    gap = arb("1.59808e-8")  # Certified for the exact gauge mean in the companion.
    times = [Fraction(t,256) for t in (1,4,8,12,15,16)]
    weights = [prod(v/(v-u) for v in times if v!=u)/u for u in times]
    for n in range(1,len(times)+1):
        assert sum(w*t**n for w,t in zip(weights,times)) == int(n==1)

    def signal(nodes, weights, calibration):
        order = len(nodes)+1
        remainder = (calibration*shell_bound * K**(order-1)
                     * (order+ball(max(nodes))*K)*ball(prod(nodes))/factorial(order))
        return (gap-remainder)*eta/(calibration*ball(sum(map(abs,weights))))

    one = Fraction(1,2**30)
    single = signal([one], [1/one], length)
    six = signal(times, weights, length)
    # Three Cartesian probabilities can also be combined in data analysis.
    cartesian = signal(times, weights, sum(map(abs,normal)))
    assert single > arb("6.8e-18") and six > arb("1.50e-11") and cartesian > arb("1.38e-11")
    print("Probability gaps per off-center fraction: one duration, six, Cartesian",
          single.str(9), six.str(9), cartesian.str(9))

    response = study["native"]["plane_response"]
    blocks = [response(acb.zeta_zero(n).imag) for n in labels]
    center = response(acb.zeta_zero(11431).imag)
    def probability(gate):
        return (event*gate*rho*gate.conjugate().transpose()).trace().real
    differences = []
    for t in times:
        U = ordered["flavor"]["identity"]
        for block,count in zip(blocks,ticks):
            U = ordered["evolution"](block,ball(t)*count/10**12,study["coefficients"])*U
        target = ordered["evolution"](center,ball(t),study["coefficients"])
        differences.append(probability(U)-probability(target))
    print("Listed integer six-segment example, full calibration box, probability differences",
          *[v.str(9) for v in differences])
    assert differences[-1] > arb("0.00007803")


if __name__ == "__main__":
    verify()

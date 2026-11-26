"""Two eliminated scalar modes suffice, and one cannot realize the alternative scale."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

root = Path(__file__).resolve().parent
clock = run_path(str(root/"ch39-coherent-clock-normalization.py"))
base, B, root_cell = [clock[k] for k in ("base","B","root_cell")]


from functools import cache


@cache
def base_weights():
    return clock["source"]["tables"]()[0]


def calibration():
    weights = base_weights()
    ids = sorted(weights)
    prefix = base["data"]("ch39-tick-separators.json")["prefix"]
    z = [clock["source"]["amb"]["upper"]["height_cell"](prefix[i-1][0]) for i in ids]
    votes = [weights[i]*x for i,x in zip(ids,z)]
    m1 = sum(votes,B(0))
    m2 = sum((v*x for v,x in zip(votes,z)),B(0))
    target = m2/(m1*m1)
    p = [v/m1 for v in votes]
    remainder = sum(p[2:],B(0))
    u, v = map(root_cell,p[:2])
    k = Q(17,20)
    # Rows of U at the two pivots are (3,0), (-1/3,3).
    # Their inner product is -1. Both have scalar product 1 with h.
    h = (Q(1,3),Q(10,27))
    assert 3*h[0]==1 and -h[0]/3+3*h[1]==1
    gamma = sum(y*y for y in h)
    assert (1+gamma)*k*k < 1
    gain = B(2)*u*v + B(2*k)*(u+v)*root_cell(remainder)
    gain -= B(gamma*k*k)*(remainder - sum((a*a for a in p[2:]),B(0))/remainder)
    upper = B(1)+gain
    assert Q("2.410")<target.lo<=target.hi<Q("2.411")
    assert Q("2.4613")<upper.lo<=upper.hi<Q("2.4614")
    theta = (target-B(1))/gain
    assert 0<theta.lo<=theta.hi<1
    # E of the microscopic block is at most 32*s^2*theta*sum(z_i |U_i|^2).
    cost = z[0]*B(9)+z[1]*B(Q(82,9))
    cost += B(gamma*k*k)*sum((zi*pi for zi,pi in zip(z[2:],p[2:])),B(0))/remainder
    assert (B(32)*theta*cost).hi<4
    # Rephase only the first height relative to a fixed coherent state.
    # Every cyclic trace of the action with the diagonal clocks stays equal.
    phase = B(1)-gain/(B(2)*u*(v+B(k)*root_cell(remainder)))
    assert Q("-.202759") < phase.lo <= phase.hi < Q("-.202758")
    shift = (B(1)-phase)*theta*root_cell(z[0]*z[1])
    assert Q(".00390758") < shift.lo <= shift.hi < Q(".00390760")
    assert Q(1,256) < shift.lo <= shift.hi < Q(1,128)
    return ids, z, p, remainder, k, gamma, theta


def main():
    calibration()
    print("A real equal-weight probe on heights 1 and 2 separates the rephased actions at 9 bits.")
    print("Target ratio between 2.410 and 2.411; one mediator gives at most 2.")
    print("Two modes reach 2.4613 before interpolation; a positive definite action realizes the target.")
    print("For every b, interior scale 2^(-b-1) keeps the microscopic commutator norm below 2^-b.")


if __name__=="__main__":
    main()

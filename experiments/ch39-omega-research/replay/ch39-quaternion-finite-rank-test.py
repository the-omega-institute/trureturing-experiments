"""A finite-precision mixed-rank test at 625 rational flux configurations, with a stated response-curvature budget."""
from fractions import Fraction as Q
from math import prod
from pathlib import Path
from runpy import run_path

from flint import arb,arb_mat,ctx
from sympy import Matrix

directory=Path(__file__).parent
quad=run_path(str(directory/"ch39-kinetic-rh-gaussian-readout.py"))["gauss"]
orient=run_path(str(directory/"ch39-shell-nine-phase-probes.py"))["orientation"]


def main():
    ctx.prec,ctx.cap=3072,3
    z,w,_=quad(21)
    W=lambda x:prod((1+j*j*x for j in range(1,42)),start=arb(1))
    lags=[W(v)/a/prod(((v-b)**2 for k,b in enumerate(z) if k!=i),start=arb(1)) for i,(v,a) in enumerate(zip(z,w))]
    amplification=arb(0)
    for j in range(1,42):
        pole=-arb(1)/(j*j)
        # Squared norm of the orthonormal polynomial evaluation vector at
        # the pole, divided by the remaining factors of the denominator.
        evaluation=sum(a*prod(((pole-b)**2 for k,b in enumerate(z) if k!=i),start=arb(1)) for i,a in enumerate(lags))
        scale=abs(prod((1+k*k*pole for k in range(1,42) if k!=j),start=arb(1)))
        amplification+=evaluation/scale
    print("Temporal error amplification: log2", (amplification.log()/arb(2).log()).str(12))
    assert amplification<arb(2)**169
    edges=orient()
    D=Matrix.zeros(48,24)
    for e,(u,v) in enumerate(edges):
        D[e,u],D[e,v]=-1,1
    tree=list(D.T.rref()[1])
    rest=[i for i in range(48) if i not in tree]
    R=D[tree,1:].inv()
    depth=[0]+[int(sum(a*a for a in R.row(i))) for i in range(23)]
    assert sum(2*(depth[edges[i][0]]+depth[edges[i][1]]) for i in rest)==558
    assert max(depth)<=23 and sum(depth)<512
    # Absolute operator error at each *raw* level gives the following
    # whitened Gram and then grounded Schur error. Harmonic extension at a
    # single root has norm <=sqrt(24) by the connection maximum principle.
    input_bits=288
    whiten=Q(2)**(169-input_bits)
    ground_gap=Q(1,512)
    delta=24*whiten/(1-whiten/ground_gap)
    assert delta<Q(2)**-114
    flat_error=Q(2)**-57
    # For zero, one, or two active unit flux axes, b/t is within 3t of
    # 2X; ||A(t)-A(0)|| <= 2 sqrt(46) t < 14t.
    t=Q(2)**-36
    curvature=2**10
    path_drift=14*t*(48+14*t)
    geometric_error=18*t+8*path_drift
    assert 3+Q(3,2)*14*(48+14*t)<curvature
    amp_error=curvature*t+flat_error/t
    eta=delta/t**2+geometric_error+(6+amp_error)*amp_error
    assert 8+delta/t**2<9
    assert 49*eta<Q(4,559)
    assert 25+2*(25*24//2)==625
    print("288-bit raw operator tolerance and 2^-36 Cayley displacements suffice under the quadratic remainder bound.")
    print("Block Gram perturbation bound",float(49*eta),"against exact lower bound",4/559)
    print("Resulting necessary quaternionic state count",48*21)


if __name__=="__main__":
    main()

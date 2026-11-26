"""Zero selection with zero diagonal and permutation-invariant positive edge magnitudes."""
from fractions import Fraction as F
from itertools import permutations
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_series, ctx

shell = run_path(str(Path(__file__).with_name("ch39-refined-mass-obstruction.py")))
ball = shell["ball"]
numerators = (2036572,2036508,2036409,2036367,2036277,2036235)


def score(tau, h):
    s = tau*tau
    u = 2*(1-s)**2/(1+s)**2
    v = 4*tau*(1-s)/(1+s)**2
    k = h*h
    q0 = 2*u*(1+k)+k*k
    q1 = 2*u+k*k+k*(-u+arb(3).sqrt()*v)
    return (2+k+q0.sqrt()).sqrt()+(2+k+q1.sqrt()).sqrt()


def verify():
    ctx.prec = 224
    nodes, edges = shell["geometry"]()
    coords = ((21,51,69), (33,33,75), (39,39,69), (19,61,61))
    classes = [set() for _ in range(4)]
    for perm in permutations(range(3)):
        a,b,v,t = [nodes.index(tuple(p[i] for i in perm)) for p in coords]
        for group, j in zip(classes, (b,v,t)):
            group.add(tuple(sorted((a,j))))
    pole = nodes.index((51,51,51))
    for i,j,_,k in edges:
        if pole in (i,j) and shell["LOCAL"]["squared"](nodes[i], nodes[j]) == 204:
            classes[3].add((i,j))
    assert all(len(group) == 6 for group in classes)
    known = set.union(*classes)
    assert len(known) == 24 and known <= {(i,j) for i,j,_,k in edges}
    marked = sorted({v for e in set.union(*classes[:3]) for v in e})
    assert len(marked) == 15 and all(i not in marked and j not in marked for i,j in classes[3])
    charges = {(i,j) for i,j,_,k in edges if k}
    assert charges == classes[0]
    # Check the Bloch expression against the eigenvalues of the marked block.
    h = F(numerators[0],10**8)
    gains = (F(1),F(1),h,F(2))
    weighted = [(i,j,g,k) for i,j,_,k in edges for group,g in zip(classes,gains) if (i,j) in group]
    tau = arb(1)/64
    phase = (1+acb(0,tau))/(1-acb(0,tau))
    mat = shell["matrices"](nodes,weighted,phase)
    block = acb_mat([[mat[i,j] for j in marked] for i in marked])
    eigen = sorted((v.real for v in block.eig(multiple=True)),key=lambda v:float(v.mid()),reverse=True)
    assert eigen[0] > eigen[1] > eigen[2]
    assert (eigen[0]+eigen[1]).overlaps(score(tau,ball(h)))

    tau_box = arb(0).union(arb(1)/28)
    h_box = arb(".020").union(arb(".021"))
    derivative = 2*score(arb_series([tau_box,1],prec=3), h_box)[2]
    assert -6 < derivative < -2
    degree = [0]*len(nodes)
    for i,j,_,_ in edges:
        if (i,j) not in known:
            degree[i] += 1
            degree[j] += 1
    cost = 6*max(degree)*ball(F(1,10**22))
    zeros = [1/(2*acb.zeta_zero(n).imag) for n in range(11484,11492)]
    assert all(a>b for a,b in zip(zeros,zeros[1:]))
    for i,num in enumerate(numerators,1):
        h = ball(F(num,10**8))
        assert h_box.contains(h)
        left, peak, right = [score(t,h) for t in zeros[i-1:i+2]]
        margin = min(peak-left,peak-right)-cost
        assert margin > arb("6.68e-18")
        print(11484+i, "symmetric link numerator", num, "global margin >", margin.lower().str(10))


if __name__ == "__main__":
    verify()

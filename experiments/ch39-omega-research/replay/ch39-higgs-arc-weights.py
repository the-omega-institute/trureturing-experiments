"""Small-arc moments with the three weights fixed by minimum Higgs transport."""
from fractions import Fraction as F
from itertools import permutations
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx

sibling = Path(__file__).with_name
opt = run_path(str(sibling("ch39-higgs-exact-optimum.py")))
arc = run_path(str(sibling("ch39-arc-step-count.py")))
q = arc["q"]


def mode_weights():
    ctx.prec = 420
    data = opt["pencils"]()
    original = [q(F(n,10**12)) for n in (190228603135,393680674736,194536323233,
                                         508983786061,701252351895,57318735330)]
    state = [x.mid() for x in original]
    for _ in range(5):
        f, j = opt["equations"](state, data)
        state = [(x-d).mid() for x,d in zip(state, j.solve(f).entries())]
    radius = q("1e-60")
    box = [x+arb(0,radius) for x in state]
    assert all(abs(v-c) < arb(2)**-33 for v,c in zip(box,original))
    f, j = opt["equations"](state, data)
    inv = arb_mat([[z.mid() for z in row] for row in j.inv().tolist()])
    _, whole = opt["equations"](box, data)
    identity = arb_mat([[int(i==j) for j in range(6)] for i in range(6)])
    slope = identity-inv*whole
    defect = inv*f
    assert all(abs(defect[i,0])+radius*sum(abs(slope[i,j]) for j in range(6)) < radius
               for i in range(6))
    return [box[4],box[5],1-box[4]-box[5]]


def verify():
    weights = mode_weights()
    assert min(weights)/3 > q("0.01831")
    reps = ((17,17,17),(11,11,25),(13,13,23),(1,5,29),(7,17,23))
    sizes = [len(set(permutations(v))) for v in reps]
    assert sizes == [1,3,3,6,6]
    for v in reps:
        assert sum(j*j for j in v)==867
    assignments = (
        ([1]*4, ((4,2,2,4),(0,3,4,2),(2,4,4,1)), (501,771,838,882,902)),
        (sizes, ((2,2,4,4,2),(4,1,0,4,4),(4,4,3,3,4)), (495,564,823,862,901)),
    )
    for counts,table,seeds in assignments:
        mixed = [arb(0) for _ in range(5)]
        for weight,row in zip(weights,table):
            for n,group in zip(counts,row):
                mixed[group] += weight*q(F(n,sum(counts)))
        assert all(x > 0 for x in mixed)
        for zero,(bounds,_) in arc["dyadic"]["DATA"].items():
            targets = [arb(q(F(2*b+1,2*10**20)),q(F(1,2*10**20))) for b in bounds]
            initial = [q(F(n,1000)) for n in seeds]
            arc["upper_certificate"](initial,targets,2**53,mixed)
        print("All six zero heights and five target harmonics verified for labels:",sum(counts))


if __name__ == "__main__":
    verify()

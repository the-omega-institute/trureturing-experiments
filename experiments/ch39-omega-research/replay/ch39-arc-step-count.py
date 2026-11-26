"""Two-sided finite lengths for independent reversible arc steps with a shared state."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx

path = Path(__file__).with_name
fifth = run_path(str(path("ch39-fifth-shell-harmonic.py")))
dyadic = run_path(str(path("ch39-dyadic-color.py")))
q = fifth["q"]
ctx.prec = 420


def equations(state,target,functions,fixed=()):
    weights, nodes = (fixed,state) if fixed else (state[:3],state[3:])
    values,jac = [],[]
    for i,goal in enumerate(target):
        terms = [functions(i+int(bool(fixed)),v) for v in nodes]
        values.append([sum(w*term[0] for w,term in zip(weights,terms))-goal])
        columns = [] if fixed else [term[0] for term in terms]
        jac.append(columns+[w*term[1] for w,term in zip(weights,terms)])
    return arb_mat(values),arb_mat(jac)


def refine(state,target,functions,fixed=()):
    for _ in range(10):
        res,jac = equations(state,target,functions,fixed)
        delta = jac.solve(res)
        state = [(v-d).mid() for v,d in zip(state,delta.entries())]
    return state


def upper_certificate(state,targets,N,fixed=()):
    x = q(F(1,2**27))
    phase = acb((1-x*x)/(1+x*x),2*x/(1+x*x))
    d = 1-phase.real
    ratios = [(1-(phase**j).real)/d for j in range(max(6, len(targets)+int(bool(fixed))))]
    def values(j,v):
        atom = v**q(F(1,N))
        base = 1-ratios[j]*(1-atom)
        value = base**N
        return value,value*ratios[j]*atom/(v*base)
    center = refine(state,[z.mid() for z in targets],values,fixed)
    rad = q("1e-12")
    box = [s+arb(0,rad) for s in center]
    weights,nodes = (fixed,box) if fixed else (box[:3],box[3:])
    assert all(w > 0 for w in weights)
    assert (1-d)**N < nodes[0] and all(a < b for a,b in zip(nodes,nodes[1:])) and nodes[-1] < 1
    res,jac = equations(center,targets,values,fixed)
    size = len(center)
    inv = arb_mat(size,size,[v.mid() for v in jac.inv().entries()])
    assert not inv.det().contains(0)
    _,whole = equations(box,targets,values,fixed)
    slope = arb_mat([[int(i==j) for j in range(size)] for i in range(size)])-inv*whole
    defect = inv*res
    assert all(abs(defect[i,0])+rad*sum(abs(slope[i,j]) for j in range(size)) < rad for i in range(size))
    return nodes


def verify():
    def ideal(i,v):
        k=i*i
        return v**k, k*v**(k-1) if k else arb(0)
    for zero,(bounds,seeds) in dyadic["DATA"].items():
        targets = [arb(1)]+[arb(q(F(2*b+1,2*10**20)),q(F(1,2*10**20))) for b in bounds]
        start = [q(v) for v in (".018",".44",".542",".525",".827",".900")]
        state = refine(start,[z.mid() for z in targets],ideal)
        v,b,c = state[3:]
        base = 1-q(F(1,2**53))
        nominal = int((v.log()/base.log()).floor().unique_fmpz())
        lower,upper = nominal-300000,nominal+25000
        assert 2**52 < lower < upper < 2**53
        cut = (base**lower).mid()
        # Pick a dyadic value provably below the single-step product lower bound.
        cut -= q("1e-28")
        assert cut < base**lower and cut < b < c < 1
        constraints = [(cut,False),(b,False),(b,True),(c,False),(c,True)]
        poly = fifth["constrained_polynomial"](constraints,fifth["POWERS"],-1)
        moment_box = fifth["moment_boxes"](targets[1:],fifth["POWERS"])
        margin = sum(a*m for a,m in zip(poly,moment_box))
        assert margin > 0
        upper_certificate(state,targets,upper)
        for counts,row in zip(dyadic["COUNTS"],seeds):
            old_d = q(dyadic["decrement"](1))
            initial = [(1-old_d*q(F(k,2**52)))**(2**56) for k in row]
            weights = [q(F(c,sum(counts))) for c in counts]
            upper_certificate(initial,targets[1:],2**53,weights)
        print(zero,lower,upper,"steps; margin",margin.str(6))


if __name__ == "__main__":
    verify()

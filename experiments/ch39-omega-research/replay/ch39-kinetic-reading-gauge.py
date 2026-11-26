"""Integer replay for the third-pole transition over the gauge rounding box."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

root=Path(__file__).resolve().parent
cross=run_path(str(root/"ch39-kinetic-reading-crossing.py"))
B=cross["B"]
rounded=cross["rounded"]


def pi_cell():
    def atan_inverse(v):
        t=sum(Q((-1)**k,(2*k+1)*v**(2*k+1)) for k in range(48))
        return B(t,t+Q(1,97*v**97))
    return rounded(B(16)*atan_inverse(5)-B(4)*atan_inverse(239))


def exp_cell(x):
    def at(t):
        assert abs(t)<2
        total=term=Q(1)
        for k in range(1,65):
            term=term*t/k
            total+=term
        err=2*abs(term*t)/65
        return rounded(B(total-err,total+err))
    return B(at(x.lo).lo,at(x.hi).hi)


def input_boxes(model,extra_bits=0):
    pi2=pi_cell()**2
    gs=[Q(".467453"),Q(".63811"),Q("1.0583")]
    half=[Q(".0000005"),Q(".000005"),Q(".00005")]
    half=[v/2**extra_bits for v in half]
    def delta(row,sign):
        diffs=[8*(1/(g+sign*(1 if a<0 else -1)*h)**2-1/g**2) for g,h,a in zip(gs,half,row)]
        return rounded(pi2*B(sum(a*d for a,d in zip(row,diffs))))
    k=(Q(5,22),Q(-3,22),Q(-2,22))
    f=(Q(-115,2),Q(333,2),Q(-109))
    # Independent enclosure of q2 and the gauge target also includes every correlated corner.
    qbounds=[rounded(model.q2*exp_cell(delta(k,i))) for i in (-1,1)]
    targets=[rounded(model.target*exp_cell(delta(f,i))) for i in (-1,1)]
    return B(qbounds[0].lo,qbounds[1].hi),targets


WINDOWS=((".00080758359582",".00080758359587"),(".00084565135521",".00084565135527"))


def main():
    model=cross["ReadingBounds"]()
    model.q2,targets=input_boxes(model)
    records=cross["base"]["data"]("ch39-kinetic-reading-gauge-roots.json")
    windows=[]
    for target,samples,(x,y) in zip(targets,records,WINDOWS,strict=True):
        model.target=target
        a,b=Q(x),Q(y)
        assert [(Q(row["c"]),row["mass"],row["sign"]) for row in samples]==[(a,None,-1),(b,None,1)]
        for row in samples:model.verify_sample(row)
        windows.append(B(model.reading(b).lo,model.reading(a).hi))
    assert Q("1.9206139144")<windows[0].lo<windows[0].hi<Q("1.9206139147")
    assert Q("1.8371077014")<windows[1].lo<windows[1].hi<Q("1.8371077017")
    print("Critical kinetic reading over gauge rounding: (1.8371077014, 1.9206139147).")
    print("Each extreme is approached to better than 0.0000000003 by an actual input corner.")


if __name__=="__main__":
    main()

"""给定第二处族流增量时第一处的锐上界及所有相容编号。"""
from pathlib import Path
from runpy import run_path

from flint import arb

four = run_path(str(Path(__file__).with_name("ch39-four-current-probes.py")))
poles,mat,q,form,j = (four[key] for key in ("poles","mat","q","form","family"))


def check(box,geometry,h,k0,cap,feasible):
    # Exclude even under measurement tolerance; prove existence for exact data.
    records = four["exchanges"](box,geometry,h,k0,q(0) if feasible else q("2e-8"))
    if records is None:
        return False
    (_,g1,g2),(low,high) = records
    K = g2-g1*cap
    if not (high>0 and low>0 and mat["positive"](K)):
        return False
    inv = K.inv()
    u,v = g1*j,g2*j
    alpha = form(v,inv,u)/(form(v,inv,v)+high)
    a = alpha.mid()
    r = u-v*a
    val = form(r,inv,r)+a**2*high
    if not feasible:
        return cap*val < low
    correction = (form(v,inv,r)-a*high)**2/(form(v,inv,v)+high)
    independent = u[0,0]*v[1,0]-u[1,0]*v[0,0]
    return cap*(val-correction)>low and not alpha.contains(0) and not independent.contains(0)


def cover(lo,hi,geometry,h,k0,cap,yes):
    pending = [(lo,hi)]
    count = 0
    while pending:
        a,b = pending.pop()
        try:
            ok = check(a.union(b),geometry,h,k0,cap,yes)
        except ZeroDivisionError:
            ok = False
        if ok:
            count += 1
        else:
            if b-a < q("1e-10"):
                raise ArithmeticError(f"undecided box {a} to {b}, expected {yes}")
            mid = ((a+b)/2).mid()
            pending.extend(((a,mid),(mid,b)))
    return count


def verify():
    time = poles["study"]["nominal"][2]
    for case,(logs,first) in enumerate(zip(poles["retune"]["log_sets"],(10946,10946,10949),strict=True)):
        labels = (1,first-1,first,11430,11432)
        targets = [poles["study"]["comp"]["model"]["target"](n,time) for n in labels]
        geometry,h,k0,ks,cap = four["setup"](logs,time,targets)
        counts = [cover(ks[0].lower(),ks[1].upper(),geometry,h,k0,cap,False),
                  cover(ks[2].lower(),ks[3].upper(),geometry,h,k0,cap,True)]
        p,a,_,_,_,refj,_ = geometry[2]
        t32 = p*ks[4]+mat["hypercharge"]*h
        g32 = (four["eye"]-t32*a).inv()*t32
        assert g32[2,2]/ks[4] > refj*(1+q("2e-8"))
        print("case",case,"first admissible label",first,"last",11431,
              "count",11432-first,"continuous boxes",counts,flush=True)


if __name__=="__main__":
    verify()

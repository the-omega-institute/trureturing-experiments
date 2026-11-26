"""Upper supporting line of the squared-height penalty; certify the later competitor by arc intervals."""
from pathlib import Path
from flint import arb,acb
import runpy

d=runpy.run_path(str(Path(__file__).with_name("ch39-shell-earlier-zero-dominance.py")))
observe=d["obs"]
T=acb.zeta_zero(719).imag
nums=(1,0,1,3,2,2)
values=[arb(v)/3 for v in nums]
def energy(h):
    return sum((x-y)**2 for x,y in zip(observe(h)[1],values))
u=energy(T)
best=arb(0)
found={}
def candidate(idx):
    global best
    if idx in found:
        return
    t=acb.zeta_zero(idx).imag
    b=(u-energy(t))/(t*t-T*T)
    found[idx]=b
    if b>best:
        best=b.lower()
def scan(a,b):
    box=a.union(b)
    # Discard only when the incumbent beats the whole continuous interval.
    if energy(box)>u-best*(a*a-T*T):
        return
    if b-a>arb(1)/8:
        mid=(a+b)/2
        scan(a,mid)
        scan(mid,b)
        return
    n0=a.zeta_nzeros().unique_fmpz()
    n1=b.zeta_nzeros().unique_fmpz()
    assert n0 is not None and n1 is not None
    for idx in range(int(n0)+1,int(n1)+1):
        if idx>719:
            candidate(idx)
def run():
    global best
    a=int(T.floor().unique_fmpz())
    # Start above the target and separately certify that this first integral endpoint has the same count.
    a+=1
    assert arb(a).zeta_nzeros().unique_fmpz()==719
    while best<=0 or best*(a*a-T*T)<=u:
        scan(arb(a),arb(a+1))
        a+=1
        if a>5000000:
            raise ArithmeticError("search cap without a tail proof")
    winner=max(found,key=lambda n:float(found[n].mid()))
    assert found[winner]>0 and all(found[winner]>z for n,z in found.items() if n!=winner)
    print("sharp lower height penalty",winner,found[winner].str(16),"tail begins",a,"zero evaluations",len(found))
    assert winner==10044 and arb(1)/51**6<found[winner]<arb(1)/51**5
    return winner,found[winner]
if __name__=="__main__":
    run()

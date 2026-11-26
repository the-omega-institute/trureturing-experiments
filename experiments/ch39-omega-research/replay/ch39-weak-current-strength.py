"""Enclose the normalized strength of the lowest weak-mode pair's harmonic current."""
from itertools import product
from pathlib import Path
import runpy
from flint import acb, acb_mat, arb

doublet = runpy.run_path(str(Path(__file__).with_name("ch39-weak-doublet-current.py")))
weak = doublet["weak"]
q = weak["q"]
reps = [o[0] for o in doublet["orbits"] if len(o)==3]
omega = acb(-q(1)/2,q(3).sqrt()/2)
assignment = {}
for a,i in enumerate(reps):
    for phase in (acb(1),omega,omega.conjugate()):
        assignment[i]=a,phase
        i=doublet["perm"][i]

short=weak["short"]
sh=weak["old"]["sh"]
raw=[(i,j,short["r"]*(sh["rho"]*sh["theta"](d)).sqrt(),short["band"](d),d)
     for i,j,d,_ in sh["edges"]]

def assess(x, factors):
    H=acb_mat(6,6)
    data=[]
    for i,j,base,band,d in raw:
        L=base*factors[band]
        w=(L/2)**2
        y=x*L/2
        sn,cs=y.sin_cos()
        off=y/sn
        for v in (i,j):
            if v in assignment:
                a,_=assignment[v]
                H[a,a]+=off*cs
        if i in assignment and j in assignment:
            a,p=assignment[i]
            b,t=assignment[j]
            term=off*p.conjugate()*t
            H[a,b]-=term
            H[b,a]-=term.conjugate()
        # Integrals of sin(y s)^2/sin(y)^2 and of the product with
        # the solution reflected about s=1/2 on a unit edge.
        diag=(1-sn*cs/y)/(2*sn**2)
        cross=(sn/y-cs)/(2*sn**2)
        data.append((i,j,w,off,diag,cross,d))
    # One circular eigenvector per nontrivial cyclic character. Normalization
    # at the first orbit is legal whenever the complementary minor is invertible.
    small=acb_mat([[H[a,b] for b in range(1,6)] for a in range(1,6)])
    rest=small.solve(acb_mat([[-H[a,0]] for a in range(1,6)]))
    z=[acb(1)]+[rest[i,0] for i in range(5)]
    vals = [acb(0) for _ in range(19)]
    for i,(a,phase) in assignment.items():
        vals[i] = z[a]*phase
    gram_lo=gram_hi=energy_lo=energy_hi=volume=q(0)
    clocks = {724:[q(0),q(0)],804:[q(0),q(0)]}
    for i,j,w,off,diag,cross,d in data:
        product=vals[i].conjugate()*vals[j]
        square_sum=abs(vals[i])**2+abs(vals[j])**2
        value=w*(diag*square_sum+2*cross*product.real)
        # On symmetric and antisymmetric endpoint data the interpolant has
        # squared mass at least 1/6 of the sum of the two endpoint squares.
        floor=w*(vals[i].abs_lower()**2+vals[j].abs_lower()**2)/6
        gram_lo+=max(floor.lower(),value.lower())
        gram_hi+=value.upper()
        flow=product.imag
        if i in assignment and j in assignment:
            a,u=assignment[i]
            b,v=assignment[j]
            if a==b:
                flow=(u.conjugate()*v).imag*abs(z[a])**2
        e=off*flow
        energy_lo+=e.abs_lower()**2
        energy_hi+=e.abs_upper()**2
        if d in clocks:
            clocks[d][0]+=e.abs_lower()**2
            clocks[d][1]+=e.abs_upper()**2
        volume+=w
    assert gram_lo>0 and energy_lo>0
    total=(volume*4*energy_lo/gram_hi**2).lower().union(
           (volume*4*energy_hi/gram_lo**2).upper())
    return [total]+[(lo/energy_hi).lower().union((hi/energy_lo).upper()) for lo,hi in clocks.values()]

def check():
    lo,hi=q("1.19576"),q("1.21744")
    def split(box,n):
        a,b=box.lower(),box.upper()
        return [a+(b-a)*q(k).union(q(k+1))/n for k in range(n)]
    vals=[]
    for a,b,c in product(*(split(f,4) for f in short["fsbox"][:3])):
        factors=short["lengths"]([a,b,c])
        if not factors[3].overlaps(short["fsbox"][3]):
            continue
        factors[3]=factors[3].intersection(short["fsbox"][3])
        for x in split(lo.union(hi),8):
            vals.append(assess(x,factors))
    for label,col in zip(("volume times current square","724 share","804 share"),zip(*vals)):
        bound=min(v.lower() for v in col).union(max(v.upper() for v in col))
        print(label,bound.lower().str(13),bound.upper().str(13))
if __name__=="__main__":
    check()

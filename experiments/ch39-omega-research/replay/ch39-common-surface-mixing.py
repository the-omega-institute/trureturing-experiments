"""Local mixing on the flag manifold with fixed quark ratios and fixed original curvature."""
from itertools import permutations, product
from pathlib import Path
import runpy
from flint import arb,acb

common=runpy.run_path(str(Path(__file__).with_name("ch39-common-two-step-spectra.py")))
fit=common["fit"]
pairs=((0,1),(1,2),(2,0))
cap=arb("1.77")


def certify():
    walls=common["shared_surface"]()
    spectra=[fit["symbols"](-fit["block"](k,fit["ii"],fit["oo"])) for k in walls]
    bonds=[next(z for i,j,d,z in fit["s"]["edges"] if (i,j,d)==(6,fit["outer"][end],292))
           for end in (1,2)]
    bare=next(z for i,j,d,z in fit["s"]["edges"] if (i,j,d)==(11,12,392))
    i,j=fit["W"].index(11),fit["W"].index(12)
    y=-walls[0][i,j]
    assert (bare*y.conjugate()).imag>0 and arb(".034")<(y+arb(".01")*bare).arg()<arb(".038")
    # Ignore the common positive sqrt(b_i b_j) in each row and the orientation sign
    # between sectors: neither changes the rank of the three complex derivatives.
    terms=[]
    for sector,end in ((0,1),(0,2),(1,1)):
        t=spectra[sector]
        d=[v*v for v in t]
        z=bonds[end-1]
        col=[]
        for h,l in pairs:
            gap=abs(d[l])**2-abs(d[h])**2
            coeff=z*fit["omega"]**(-end*l)*(t[h]+t[l])*d[l].conjugate()/gap
            conj=z.conjugate()*fit["omega"]**(end*h)*(t[h]+t[l]).conjugate()*d[h]/gap
            col.append((coeff,conj))
        terms.append(col)
    # det is a cubic homogeneous polynomial; collect terms before any box evaluation.
    powers={}
    for perm in permutations(range(3)):
        sign=(-1)**sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3))
        for choices in product((0,1),repeat=3):
            degree=[0,0,0]
            val=acb(sign)
            for row,(h,l) in enumerate(pairs):
                branch=choices[row]
                degree[(l,h)[branch]]+=1
                val*=terms[perm[row]][row][branch]
            key=tuple(degree)
            powers[key]=powers.get(key,acb(0))+val
    def bound(box,depth):
        b=[arb(1)]+[a.union(c) for a,c in box]
        val=sum(z*b[1]**j*b[2]**k for (i,j,k),z in powers.items())
        if abs(val)>arb(".030"):
            return 1
        assert depth<24
        j=depth%2
        a,c=box[j]
        mid=(a+c)/2
        sub=[]
        for part in ((a,mid),(mid,c)):
            new=list(box)
            new[j]=part
            sub.append(bound(new,depth+1))
        return sum(sub)
    print("scaled complex minor lower bound 0.030; ratio boxes",bound([(1/cap,cap)]*2,0))


if __name__=="__main__":
    certify()

"""Mass-preserving CP breaking on one 292 diagonal of the common degree-three surface."""
from itertools import product
from pathlib import Path
import runpy
from flint import arb,acb

common=runpy.run_path(str(Path(__file__).with_name("ch39-common-two-step-spectra.py")))
fit=common["fit"]
pairs=((0,1),(1,2),(2,0))
cap=arb("1.77")
assert (arb(793)/473)*arb(".04822").exp()<cap


def coefficient(end,e,bond):
    two=[z*z for z in e]
    gap=abs((abs(two[0])**2-abs(two[1])**2)*(abs(two[1])**2-abs(two[2])**2)*
            (abs(two[2])**2-abs(two[0])**2))
    # Eigenvectors of the circulant in fit.symbols have components omega^(-i*j).
    off=[[bond*(e[i]+e[j])*fit["omega"]**(-end*j)/3 for j in range(3)] for i in range(3)]
    def calc(r1,r2):
        b=[arb(1),r1,r2]
        z=acb(1)
        for i,j in pairs:
            z*=b[j]*off[i][j]*two[j].conjugate()+b[i]*off[j][i].conjugate()*two[i]
        return z.imag/(r1*r2*gap)
    low,high=(arb(".0016"),arb(".0072")) if end==1 else (arb("-.0072"),arb("-.0016"))
    limits=[(1/cap,cap)]*2
    cells=0
    def bound(box,depth):
        nonlocal cells
        v=calc(*(a.union(b) for a,b in box))
        if low<v<high:
            cells+=1
            return
        assert depth<30
        j=depth%2
        a,b=box[j]
        mid=(a+b)/2
        for part in ((a,mid),(mid,b)):
            new=list(box)
            new[j]=part
            bound(new,depth+1)
    bound(limits,0)
    print("292 diagonal endpoint, signed leading coefficient interval, boxes",end,low,high,cells)


if __name__=="__main__":
    _,d=common["shared_surface"]()
    t=-fit["block"](d,fit["ii"],fit["oo"])
    e=fit["symbols"](t)
    for i,j in product(range(3),repeat=2):
        v=sum(t[0,n]*fit["omega"]**(-j*n) for n in range(3))*fit["omega"]**(-j*i)
        assert (v-sum(t[i,n]*fit["omega"]**(-j*n) for n in range(3))).contains(0)
    for end in (1,2):
        a,b=6,fit["outer"][end]
        assert not any((a,b) in tuple(zip(f,f[1:]+f[:1])) or
                       (b,a) in tuple(zip(f,f[1:]+f[:1])) for f in fit["src"]["faces"])
        bond=next(z for i,j,d,z in fit["s"]["edges"] if (i,j,d)==(a,b,292))
        coefficient(end,e,bond)

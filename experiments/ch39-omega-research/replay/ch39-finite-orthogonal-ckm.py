"""Seven original faces and a finite CKM fit at zero 719 in orthonormal ring labels."""
from pathlib import Path
import runpy
from flint import acb, arb, acb_mat

common=runpy.run_path(str(Path(__file__).with_name("ch39-common-two-step-spectra.py")))
fit,flux=common["fit"],common["flux"]
wall=fit["W"]
weights=(
    (".300470260509","1.69686796692","1.81671473275",".000001269530",
     ".311143930602",".000001014555",".105528633408",".562418160757",".074617230877",".294194577122"),
    (".949403886184",".477262503187","2.39444845504",".083830714857",
     ".385439236337",".000001023568",".671673882377",".377119355527",".069407511409",".009351266279"))
identity=acb_mat([[int(i==j) for j in range(3)] for i in range(3)])


def bond(i,j):
    a,b=sorted((i,j))
    z=next(z for u,v,_,z in fit["s"]["edges"] if (u,v)==(a,b))
    return z if i==a else z.conjugate()


def add_edge(k,i,j,w):
    z=bond(i,j)
    a,b=wall.index(i),wall.index(j)
    k[a,a]+=w
    k[b,b]+=w
    k[a,b]-=w*z
    k[b,a]-=w*z.conjugate()


def fix_phase(k,cycle,angle):
    i,j=cycle[:2]
    y=-k[wall.index(i),wall.index(j)]
    rot=acb(0,-angle).exp()
    w=-(y*rot).imag/(bond(i,j)*rot).imag
    assert w>0 and ((y+w*bond(i,j))*rot).real>0
    for a,b in zip(cycle,cycle[1:]+cycle[:1]):
        add_edge(k,a,b,w)
    return w


def transfer(model,weight):
    w=list(map(arb,weight))
    k=fit["evaluate"](model,w[:4],(0,0))[2]
    for end in (1,2):
        for j in range(3):
            add_edge(k,fit["inner"][j],fit["outer"][(j+end)%3],w[4+3*(end-1)+j])
    inn=fix_phase(k,(7,6,8),arb.pi()-arb(".006"))
    base=fix_phase(k,(11,12,18),arb(".03606860"))
    degrees,margin=flux["degrees"](k)
    assert degrees[0]==3 and degrees[1]!=[3] and margin>0
    print("200/392 stiffnesses, seven faces and eight splits",inn.str(9),base.str(9),degrees)
    t=-fit["block"](k,fit["ii"],fit["oo"])
    return t*t


def projectors(y):
    h=y*y.conjugate().transpose()
    eig=sorted((v.real for v in h.eig(algorithm="rump")),key=lambda v:float(v.mid()))
    assert 0<eig[0]<eig[1]<eig[2]
    projs=[]
    for j,lam in enumerate(eig):
        p=identity
        for i,other in enumerate(eig):
            if i!=j:
                p=p*(h-other*identity)/(lam-other)
        projs.append(p)
    return [(v/eig[2]).sqrt() for v in eig[:2]],projs


def observables(yu,yd):
    mu,p=projectors(yu)
    md,q=projectors(yd)
    def trace(m):
        return sum(m[i,i] for i in range(3))
    def overlap(i,j):
        return trace(p[i]*q[j]).real.sqrt()
    s13=overlap(0,2)
    c13=(1-s13*s13).sqrt()
    s12,s23=overlap(0,1)/c13,overlap(1,2)/c13
    quartet=trace(p[0]*q[1]*p[1]*q[2])
    delta=(quartet/(s12*s23*s13*c13**2)+s12*s23*s13).arg()
    return mu+md+[s12,s23,s13,delta,quartet.imag]


if __name__=="__main__":
    model=fit["model"](*common["angles"])
    vals=observables(*(transfer(model,w) for w in weights))
    centers=common["up"]+list(reversed(common["down"]))+list(map(arb,(".22501",".04183",".003732","1.147")))
    for value,ref in zip(vals,centers):
        assert abs(value/ref-1)<arb("1e-6")
    assert arb(".0000300")<vals[-1]<arb(".0000325")
    print("four mass ratios, s12/s23/s13, phase, J",*(v.str(10) for v in vals))

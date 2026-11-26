"""Test the lepton angles for all relative coordinate permutations and independent orbit weights."""
from functools import lru_cache
from itertools import permutations,product
from pathlib import Path
from runpy import run_path
from flint import arb,arb_mat,acb_mat

m=run_path(str(Path(__file__).with_name("ch39-shell-four-seesaw.py")))
eye=m["eye"]
d=arb_mat([[1],[1],[1]])/arb(3).sqrt()
P=eye-d*d.transpose()
indices=(0,1,2)
perms=list(permutations(indices))


def vec(x):
    return [x[i,0] for i in indices]


def unit(x):
    return x/sum(v*v for v in vec(x)).sqrt()


def dot(x,y):
    return sum(a*b for a,b in zip(x,y))


def outside(x,lo,hi):
    return x<lo or x>hi


def check(n):
    H=sum(m["matrices"](n),arb_mat(3,3))/4
    B=P*H*P
    trace=sum(B[i,i] for i in indices)
    gap=(2*sum((B*B)[i,i] for i in indices)-trace**2).sqrt()
    beta=((trace-gap)/2,(trace+gap)/2)
    avg=trace/2
    k=(d.transpose()*H*d)[0,0]-avg
    assert k>0
    U=unit((P*beta[1]-B)*H*d)
    V=unit((B-P*beta[0])*H*d)
    g=((U.transpose()*H*d)[0,0],(V.transpose()*H*d)[0,0])
    assert min(g)>0
    tr=sum(H[i,i] for i in indices)

    @lru_cache(None)
    def endpoint(level,j):
        # Exact dyadic parameters, with the two bounded eigenvalues extended to eta=0.
        if j==0:
            return beta
        eta=arb(j)/2**level
        s=k*(1-eta)/eta
        es=sorted((z.real for z in acb_mat(H+d*d.transpose()*s).eig()),key=lambda z:float(z.mid()))
        assert es[0]<es[1]<es[2]
        return tuple(es[:2])

    @lru_cache(None)
    def box(level,j):
        eta=(arb(j)/2**level).union(arb(j+1)/2**level)
        l,r=endpoint(level,j),endpoint(level,j+1)
        a,b=(x.union(y) for x,y in zip(l,r))
        light=[(1-eta)*avg+eta*v for v in (a,b)]
        light.append(tr-sum(light))
        masses=None
        if all(v>0 for v in light):
            sq=[1/v**2 for v in light[::-1]]
            if sq[2]>sq[0]:
                masses=[(sq[j+1]-sq[j])/(sq[2]-sq[0]) for j in (0,1)]
        vectors=[]
        for i,x in enumerate((a,b)):
            a0,a1=(U,V) if i==0 else (V,U)
            dx=(beta[i]-x)/g[i]
            coeff=g[1-i]*dx/(beta[1-i]-x)
            # Orthonormal coefficient norm is >=1, including over coarse intervals.
            den=(1+coeff**2+dx**2).sqrt()
            if not den>0 or not den.is_finite():
                return masses,None
            vectors.append(vec((a0+a1*coeff-d*dx)/den))
        a0,a1=vectors
        vectors.append([a0[(i+1)%3]*a1[(i+2)%3]-a0[(i+2)%3]*a1[(i+1)%3] for i in indices])
        return masses,vectors

    masslim=list(map(arb,(".02655",".034")))
    s13=(arb(".019"),arb(".025"))
    s12=(arb(".270"),arb(".350"))
    s23=(arb(".4"),arb(".64"))
    opts=[(o,p,e,mu) for o,p in product(range(2),range(6)) for e,mu,_ in perms]
    leaves=0
    def search(nu,el,options):
        nonlocal leaves
        dm,N=box(*nu)
        _,E=box(*el)
        remain=[]
        if dm is not None:
            options=[op for op in options if not outside(dm[op[0]],*masslim)]
        if N is None or E is None:
            remain=options
        else:
            memo={}
            def prob(p,i,j):
                key=(p,i,j)
                if key not in memo:
                    memo[key]=dot(E[i],[N[j][k] for k in perms[p]])**2
                return memo[key]
            for o,p,e,mu in options:
                # Column 3 is isolated; column 2 is the higher-mass member of the solar pair.
                j3,j2=((0,1) if o==0 else (2,0))
                reactor=prob(p,e,j3)
                if outside(reactor,*s13):
                    continue
                # Also intersect the denominator with its experimental range.
                den=1-reactor.intersection(s13[0].union(s13[1]))
                if outside(prob(p,e,j2)/den,*s12) or outside(prob(p,mu,j3)/den,*s23):
                    continue
                remain.append((o,p,e,mu))
        if not remain:
            leaves+=1
            return
        # Refine the mass sector first on comparable scales.
        side=0 if nu[0]<el[0]+2 else 1
        chosen=(nu,el)[side]
        l,j=chosen
        assert l<24
        for k2 in (2*j,2*j+1):
            args=[nu,el]
            args[side]=(l+1,k2)
            search(*args,remain)
    search((0,0),(0,0),opts)
    print("zero",n,"all relative permutations, 2 orderings and 6 charged labelings excluded; boxes",leaves)


if __name__=="__main__":
    for n in (25,255,719):
        check(n)

"""Exact physical scalar masses and subgroup indices of the quartic adjoint model."""
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
from runpy import run_path
import numpy as np
from sympy import Matrix, Poly, Symbol, factor_list

ROOT = Path(__file__).parent
OLD = run_path(str(ROOT/"ch39-locking-scalar-spectrum.py"))
LIFT = run_path(str(ROOT/"ch39-locking-vector-spectrum.py"))["derivative"]
P,Q,COLOR = OLD["P"],OLD["Q"],OLD["COLOR"]


def perturbations(sign):
    for f,i,j in product(range(2),range(10),range(10)):
        x,y,k = np.zeros((10,10),int),np.zeros((10,10),int),np.zeros((5,5),int)
        (x if f==0 else y)[i,j]=1
        yield x,y,k,OLD["weight"](i,j)
    for i,j in combinations(range(5),2):
        k=np.zeros((5,5),int)
        k[i,j],k[j,i]=1,sign
        w=tuple(int(r==i)-int(r==j) for r in range(5))
        yield np.zeros((10,10),int),np.zeros((10,10),int),k,w
    if sign==1:
        for i in range(4):
            k=np.zeros((5,5),int)
            k[i,i],k[i+1,i+1]=1,-1
            yield np.zeros((10,10),int),np.zeros((10,10),int),k,(0,)*5


def lift(k):
    l=np.zeros((10,10),int)
    for ij,val in LIFT({(int(i),int(j)):k[i,j] for i,j in zip(*np.nonzero(k))}).items():
        l[ij]=val
    return l


GENERATORS = [(-2, -2, -2, 3, 3), (0, 0, 0, 1, -1), (1, -1, 0, 0, 0)]


def weight_indices(weight):
    return [F(sum(a*b for a, b in zip(weight, h))**2, 2*sum(c*c for c in h))
            for h in GENERATORS]


def mass_blocks(ratio):
    """Return the real weight blocks for an integer u/v, in units of v squared."""
    blocks_by_weight = {}
    end=lift(COLOR)-2*np.eye(10,dtype=int)
    for s in (1,-1):
        kin,jac,blocks=[],[],defaultdict(list)
        for i,(x,y,k,w) in enumerate(perturbations(s)):
            dx=s*x.T@P+P@x
            parts=[dx+s*y.T@Q+Q@y,x@Q+s*P@y.T,s*x.T@Q+P@y,np.array([dx.trace()]),
                   ratio*(COLOR@k+k@COLOR-k),ratio*x@end+P@lift(k)]
            jac.append(np.concatenate([a.ravel() for a in parts]))
            kin.append(np.concatenate([a.ravel() for a in (x,y,k)]))
            blocks[w].append(i)
        J,K=np.column_stack(jac),np.column_stack(kin)
        mass,metric=J.T@J,K.T@K
        for w,rows in blocks.items():
            others=[j for key,entries in blocks.items() if key!=w for j in entries]
            assert not mass[np.ix_(rows,others)].any() and not metric[np.ix_(rows,others)].any()
            m,g=[a[np.ix_(rows,rows)] for a in (mass,metric)]
            blocks_by_weight[s, w] = (m, g)
    return blocks_by_weight


def spectrum():
    factors = Counter()
    indices = defaultdict(lambda: [F(0)]*3)
    p = Symbol("p")
    factor_cache = {}
    for (_, weight), (mass, metric) in mass_blocks(1).items():
        key = mass.shape, mass.tobytes(), metric.tobytes()
        if key not in factor_cache:
            operator = Matrix(metric.tolist()).inv()*Matrix(mass.tolist())
            polynomial = operator.charpoly(p).as_expr()
            factor_cache[key] = [(tuple(F(v) for v in Poly(f, p).monic().all_coeffs()), n)
                                 for f, n in factor_list(polynomial)[1]]
        charges = weight_indices(weight)
        for factor, count in factor_cache[key]:
            factors[factor] += count
            for i, charge in enumerate(charges):
                indices[factor][i] += count*charge
    return factors,indices


a = run_path(str(ROOT/"ch39-locking-coupling-match.py"))["a"]


def scalar_thresholds(indices):
    return [-sum(a(row[i])*a(abs(f[-1])).log() for f,row in indices.items()
                   if f[-1])/12 for i in range(3)]


def matched(scalar, epsilon, k, ratio=1):
    from flint import arb
    pi = arb.pi()
    inv = [4*pi/a(g)**2 for g in (".467453",".63811","1.0583")]
    d = [2*pi*x-s+a(35-h)/6 for x,s,h in zip(inv,scalar,(0,2,3))]
    t0 = (5*d[0]+3*d[1]-8*d[2])/67
    r0 = a("15/16")*(d[0]-d[1]-a("109/15")*t0)
    z = a("134/7")*(r0/3-epsilon.log())
    l = z.exp()/arb(2+ratio*ratio).sqrt()
    t = t0-a("63/67")*z
    g10 = l*(3+k*k).sqrt()
    mh = l*(3+k*k)/k
    vector = a("7/2")*(5*mh.log()+3*z+25*g10.log())-a("33/6")
    log_x = (8*pi**2/l**2-a("19/6")*t-6*epsilon.log()-2*pi*inv[1]+scalar[1]+vector)/6
    return l, 1000*t.exp(), log_x


def matching_checks(scalar):
    from flint import arb, acb, ctx
    ctx.prec=150
    def strength(n):
        return 1/(acb.zeta_zero(n).imag**2+arb(1)/4)
    low46, _, log_min = matched(scalar, strength(46), 1/arb(2).sqrt())
    assert low46 < a(".45") and log_min > 35
    slope = -16*arb.pi()**2/a(".45")**2 + a("19/6")*63/67 + a("42/134") + a("231/2")
    assert slope < 0
    assert matched(scalar,strength(49),arb(1))[0] > arb.pi().sqrt()
    for n in (47,48):
        l,v,lx = matched(scalar,strength(n),arb(1))
        assert arb(1)/4 < l*l/arb.pi() < 1 and lx.exp() < 1 and v*lx.exp()*strength(n)>1000
        print(n, "g_low", l.str(10), "alpha", (l*l/arb.pi()).str(10),
              "v/GeV", v.str(10), "X", lx.exp().str(10))
    shell = [p for p in product(range(30), repeat=3) if sum(c*c for c in p) == 867]
    edges = [(p, q) for p in shell for q in shell if sum((i-j)**2 for i,j in zip(p,q)) == 12]
    poles = [p for p in shell if len(set(p)) == 1]
    assert poles == [(17, 17, 17)] and poles[0] not in {p for p, _ in edges}
    assert len(shell) == 19 and len(edges) == 6 and len(set(edges)) == 6
    for eps in (arb(7)/19**4,arb(7)/(19**4-1)):
        l,v,lx = matched(scalar,eps,arb(1))
        assert arb(".662") < l < arb(".663") and arb("6.21e14") < v < arb("6.22e14")
        assert arb(".149") < lx.exp() < arb(".151")
        assert v*lx.exp()*eps > 1000 and l*l/arb.pi() < a("1/4")
        print("seven-word fit", l.str(10), v.str(10), lx.exp().str(10))
    # This interval encloses the six high-zero cubic graph readings at fixed hopping.
    assert matched(scalar,arb("4.609e-5"),arb(1))[0] > 12
    assert matched(scalar,arb(6)/19**4,arb(1))[0] > 12
    for k, bound in ((6,arb("1.5")),(8,arb(".3"))):
        l = matched(scalar,arb(k)/19**4,arb(1))[0]
        assert l > bound if k==6 else l < bound
    print("Single-zero labels 47,48 at alpha <= 1; seven four-draw events give a moderate match.")


def verify():
    counts,idx=spectrum()
    assert sum((len(f)-1)*n for f,n in counts.items()) == 424
    assert counts[(1,0)]==112 and idx[(1,0)]==[35,33,32]
    assert [sum((len(f)-1)*v[i] for f,v in idx.items()) for i in range(3)]==[125]*3
    matching_checks(scalar_thresholds(idx))
    for f in sorted(counts):
        print(tuple(map(str,f)),counts[f],tuple(map(str,idx[f])))


if __name__=="__main__":
    verify()

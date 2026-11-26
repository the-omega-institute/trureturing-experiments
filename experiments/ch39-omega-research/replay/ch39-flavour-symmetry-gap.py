"""Sharp simultaneous flavour-diagonal approximation in the two measured quark mass bases."""
from itertools import combinations
from flint import acb,arb


def ckm(s12,s23,s13,delta):
    c12,c23,c13=((1-s*s).sqrt() for s in (s12,s23,s13))
    z=acb(0,delta).exp()
    return [[c12*c13,s12*c13,s13*z.conjugate()],
            [-s12*c23-c12*s23*s13*z,c12*c23-s12*s23*s13*z,s23*c13],
            [s12*s23-c12*c23*s13*z,-c12*s23-s12*c23*s13*z,c23*c13]]


def leakage(uncertainty,phase_error):
    noise=arb(0,arb(uncertainty))
    us,cb,ub=(arb(v)+noise for v in (".22501",".04183",".003732"))
    matrix=ckm(us,cb,ub,arb("1.147")+arb(0,arb(phase_error)))
    p=[[abs(z)**2 for z in row] for row in matrix]
    weights=[sum(a*b for a,b in zip(p[i],p[j])) for i,j in combinations(range(3),2)]
    s=sum(weights)
    pairs=sum(a*b for a,b in combinations(weights,2))
    gap=3*pairs/(s+(s*s-3*pairs).sqrt())
    return (gap/(2*(1+(1-gap).sqrt()))).sqrt()


if __name__=="__main__":
    center=leakage("0","0")
    sector=leakage(".000055",".017")
    print("sharp minimax leakage, PDG centers and full RH clock square",center.str(12),sector.str(9))
    assert arb(".03612")<sector<arb(".03632")

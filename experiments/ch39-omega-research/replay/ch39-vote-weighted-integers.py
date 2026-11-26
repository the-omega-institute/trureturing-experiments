"""Integer replay of the universal conditional parent bracket on genuine zero heights."""
from fractions import Fraction as Q
from math import prod
from pathlib import Path
from runpy import run_path

root=Path(__file__).resolve().parent
base=run_path(str(root/"ch39-cutoff-thresholds.py"))
B=base["Bounds"]
A=[j*j for j in base["OFFSETS"] if j]


def round_out(x):
    size=max(abs(x.lo),abs(x.hi))
    if not size:
        return x
    step=Q(2)**(size.numerator.bit_length()-size.denominator.bit_length()-220)
    return B((x.lo//step)*step,-((-x.hi)//step)*step)


def multiply(*values):
    out=B(1)
    for v in values:
        out=round_out(out*v)
    return out


def absolute(x):
    if x.lo>0:return x
    assert x.hi<0
    return B(-x.hi,-x.lo)


def main():
    data=base["data"]
    rec=data("ch39-cutoff-threshold-inverses.json")[1]
    weights=base["exact_fit"](rec)
    prefix=data("ch39-tick-separators.json")["prefix"]
    def ordinate(i):
        h=prefix[i-1][0]
        return round_out(B(Q(2**800,(h+1)**2),Q(2**800,h*h)))
    z=[ordinate(i) for i in range(1,45)]
    assert all(u.lo>v.hi for u,v in zip(z,z[1:]))
    support=[(ordinate(i),w) for i,w in weights.items()]
    p,modes,q02,target=base["physical_inputs"]()
    zprod=multiply(*z)
    tails=[]
    for a in A:
        numer=multiply(*(v+B(Q(1,a)) for v in z))
        den=prod(abs(1-Q(k,a)) for k in A if k!=a)
        tails.append(round_out(numer/B(den)).hi)
    tol=Q(1,2**272)
    zerr=(2*zprod.hi+tol*sum(tails))*prod(A)
    assert zerr<Q("4e-55")
    mu=sum((round_out(v*w) for v,w in support),B(-zerr,zerr))
    for i,v in enumerate(z[:3]):
        other=z[:i]+z[i+1:]
        factor=multiply(*(B(1)+B(a)*v for a in A))/multiply(*(absolute(v-t) for t in other))
        tail=multiply(factor,*other).hi
        sensitivity=factor.hi*sum(h/(Q(1,a)+v.lo) for a,h in zip(A,tails))
        assert weights[i+1].lo-2*tail-tol*sensitivity>Q(".05")
    poles=[B(1)/v for v in z[:3]]
    certificates=data("ch39-vote-weighted-roots.json")
    assert [(r["mass"],r["sign"]) for r in certificates]==[(50772311235760,1),(50772311235762,-1)]
    for sample in certificates:
        mass=sample["mass"]
        q=q02/B(mass**2)
        assert q.hi<8 and (mu*poles[1]).hi<q.lo and (mu*poles[2]).lo>q.hi
        higgs=B(Q(125**2,mass**2))
        exposure=multiply(*[(q/(mu*t))**25 for t in poles[:2]])
        active=50
        for lam,n,block in zip(modes,base["MODES"],sample["roots"],strict=True):
            assert len(block)==3
            def f(y):
                point=B(y)
                totals=[round_out(w/(B(1)-point*t)) for t,w in support]
                den=multiply(*(absolute(B(1)-point*t) for t in z)).lo
                multiplier=y*prod(y+a for a in A)/den
                assert 1-y*z[-1].hi>0
                error=multiplier*(2*zprod.hi/(1-y*z[-1].hi)+tol*sum(h/(1+y/a) for h,a in zip(tails,A)))
                total=sum(totals,B(-error,error))
                return mu*point-higgs-lam+lam*(B(1)-p)*total
            for j,(l,h) in enumerate(block):
                y=B(Q(l,2**64),Q(h,2**64))
                assert y.lo>(poles[j-1].hi if j else 0) and y.hi<poles[j].lo
                assert f(y.lo).hi<0 and f(y.hi).lo>0
                mass2=mu*y
                assert mass2.lo*mass**2>1000**2
                if mass2.hi<q.lo:
                    exposure=multiply(exposure,(q/mass2)**n)
                    active+=n
                else:
                    assert mass2.lo>q.hi
        value=exposure**11
        assert active==119
        if sample["sign"]==1:assert value.lo>target.hi
        else:assert value.hi<target.lo
    print("Integer replay: all positive spectra within 2^-272 in the 44 time rows, without any height cutoff.")
    print("Vote-weighted kinetics, 119 open scalar doublets; 50772311235760 < M_F/GeV < 50772311235762.")


if __name__=="__main__":
    main()

"""Universal one-loop window for the alternative vote-weighted kinetic normalization."""
import json
import sys
from math import prod
from pathlib import Path
from runpy import run_path

from flint import arb,ctx

root=Path(__file__).resolve().parent
base=run_path(str(root/"ch39-cutoff-thresholds.py"))
ball=run_path(str(root/"ch39-fixed-poles-exposure.py"))["ball"]
A=[j*j for j in base["OFFSETS"] if j]
assert len(A)==43


def main():
    ctx.prec=512
    rec=base["data"]("ch39-cutoff-threshold-inverses.json")[1]
    weights={i:ball(v) for i,v in base["exact_fit"](rec).items()}
    prefix=base["data"]("ch39-tick-separators.json")["prefix"]
    def ordinate(i):
        h=prefix[i-1][0]
        return (arb(2)**800)/arb(h).union(arb(h+1))**2
    z=[ordinate(i) for i in range(1,45)]
    support=[(ordinate(i),w) for i,w in weights.items()]
    p,modes,q02,target=base["physical_inputs"]()
    p,q02,target=map(ball,(p,q02,target))
    modes=list(map(ball,modes))
    product_z=prod(z)
    # For interpolation on these 44 z nodes, evaluate the numerator at the
    # simple time poles -1/a to recover its partial fraction coefficients.
    # The coefficient of the constant row incurs zero error by normalization.
    tol=arb(2)**-272
    tails=[prod(1/arb(a)+v for v in z)/prod(abs(1-arb(k)/a) for k in A if k!=a) for a in A]
    error_z=(2*product_z+tol*sum(tails))*prod(A)
    assert error_z<arb("4e-55")
    mu=sum(v*w for v,w in support)+arb(0,error_z.upper())
    assert mu>0
    for i in range(3):
        # Rational cardinal function equals one at node i and zero at
        # all the other 43 anchors, and its expectation is known.
        v=z[i]
        tail=prod(1+a*v for a in A)*prod(t/abs(v-t) for k,t in enumerate(z) if k!=i)
        sensitivity=prod(1+a*v for a in A)/prod(abs(v-t) for k,t in enumerate(z) if k!=i)
        sensitivity*=sum(h/(1/arb(a)+v) for a,h in zip(A,tails))
        assert weights[i+1]-2*tail-tol*sensitivity>arb(".05")
    poles=[1/x for x in z[:3]]
    certificates=[]
    for mass,sign in ((50772311235760,1),(50772311235762,-1)):
        blocks=[]
        q=q02/mass**2
        higgs=arb(125**2)/mass**2
        assert q<8 and mu*poles[1]<q and mu*poles[2]>q
        exposure=prod((q/(mu*t))**25 for t in poles[:2])
        active=50
        for lam,n in zip(modes,base["MODES"]):
            block=[]
            def f(y):
                total=sum(w/(1-y*t) for t,w in support)
                # Remainder after interpolation of 1/(1-y*z) in the 44 calibrated rows.
                # At omitted nodes, z<=z_44, W(z)>=1, and |z-z_i|<=z_i.
                multiplier=y*prod(y+a for a in A)/prod(abs(1-y*t) for t in z)
                error=multiplier*(2*product_z/(1-y*z[-1])+tol*sum(h/(1+y/a) for a,h in zip(A,tails)))
                assert error>=0
                return mu*y-higgs-lam+lam*(1-p)*(total+arb(0,error.upper()))
            for lo,hi in zip([arb(0)]+poles,poles):
                left,right=lo.upper(),hi.lower()
                for _ in range(95):
                    mid=(left+right)/2
                    if f(mid)<0:left=mid
                    else:right=mid
                lower=left
                left,right=lo.upper(),hi.lower()
                for _ in range(95):
                    mid=(left+right)/2
                    if f(mid)>0:right=mid
                    else:left=mid
                upper=right
                assert f(lower)<0 and f(upper)>0
                block.append([int((lower*2**64).floor().unique_fmpz())-2,
                              int((upper*2**64).ceil().unique_fmpz())+2])
                y=lower.union(upper)
                assert lo<y and y<hi and mu*y*mass**2>1000**2
                if mu*y<q:
                    exposure*=(q/(mu*y))**n
                    active+=n
                else:
                    assert mu*y>q
            blocks.append(block)
        certificates.append({"mass":mass,"sign":sign,"roots":blocks})
        assert active==119
        assert sign*(exposure**11-target)>0
    receipt=root/"ch39-vote-weighted-roots.json"
    if "--write" in sys.argv:
        receipt.write_text(json.dumps(certificates)+"\n")
    else:
        assert json.loads(receipt.read_text())==certificates
    print("All finite positive realizations of the 44 time readings within 2^-272 on genuine heights:")
    print("vote-weighted kinetic mean",mu.str(18),"119 open scalar doublets.")
    print("Parent mass strictly between 50772311235760 and 50772311235762 GeV.")


if __name__=="__main__":
    main()

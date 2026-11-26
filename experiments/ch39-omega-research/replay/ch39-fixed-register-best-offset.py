"""Find the unique integer offset best separating the two dyadic scale endpoints."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path


def main():
    root=Path(__file__).resolve().parent
    lib=run_path(str(root/"ch39-fixed-register-scale.py"))
    t=lib["Table"](root)
    delta=Q(1,2**166)
    minus,errm=lib["actual_weights"](t,-delta)
    plus,errp=lib["actual_weights"](t,delta)
    assert min(minus)>errm and min(plus)>errp
    def difference(j):
        a=sum(v*u for v,u in zip(minus,t.observation(j,1-delta)))
        b=sum(v*u for v,u in zip(plus,t.observation(j,1+delta)))
        margin=t.n*(errm+errp)
        return b-a-margin,b-a+margin
    for j,bits,l,h in [(42,168,Q("1.01550e-50"),Q("1.01552e-50")),
                       (1682,24,Q(".00000021135517"),Q(".00000021135518"))]:
        lo,hi=difference(j)
        assert l<lo<hi<h
        assert 2*Q(1,2**bits)<lo<hi<2*Q(1,2**(bits-1))
        print("Offset",j,"least binary depth for disjoint error bands",bits)

    unit=1<<112
    d=[(unit*(b-a))//1 for a,b in zip(minus,plus)]
    # For normalized nonnegative weights, moving the common scale to one
    # costs at most delta/(1-delta) on each response, uniformly in j.
    pad=-((-2*delta*unit/(1-delta)-t.n*unit*(errp+errm))//1)+t.n
    def window(left,right):
        low=-pad; high=pad
        for weight,z in zip(d,t.z):
            a,b=z.numerator,z.denominator
            low_j,high_j=(right,left) if weight>=0 else (left,right)
            low+=weight*b//(b+a*low_j**2)
            high-=(-weight*b)//(b+a*high_j**2)
        return low,high
    best=1682; lower,_=window(best,best)
    cutoff=60000
    assert sum(Q(abs(n),unit)/z for n,z in zip(d,t.z))/cutoff**2+Q(pad,unit)<Q(lower,unit)
    stack=[(0,best-1),(best+1,cutoff-1)]
    while stack:
        left,right=stack.pop()
        low,high=window(left,right)
        if max(-low,high)<lower:
            continue
        assert left<right, "competing integer offset"
        mid=(left+right)//2
        stack.extend(((left,mid),(mid+1,right)))
    print("Unique largest absolute difference at j=1682, over all nonnegative integer offsets.")


if __name__=="__main__":
    main()


"""Certify the finite integer record against the eight doubled-clock zeta heights."""
import json
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

directory=Path(__file__).parent
joint=run_path(str(directory/"ch39-rh-shell-clock-step.py"))
original=run_path(str(directory/"ch39-finite-rh-pivot-mass.py"))
bits=208
den=2**bits


def tick(v):
    n=(v*den).floor().unique_fmpz()
    assert n is not None
    return int(n)


def rational_turn(angle):
    # On the better-conditioned hemisphere represent the phase as
    # sign * (1-t^2+2it)/(1+t^2). Its two coordinates are rational for dyadic t.
    x,y=angle.cos(),angle.sin()
    sign=1 if x>0 else -1
    t=sign*y/(1+sign*x)
    n=tick(t)
    u=arb(n)/den
    a=sign*(1-u*u)/(1+u*u)
    b=sign*2*u/(1+u*u)
    assert ((x-a)**2+(y-b)**2).sqrt()<arb(3)/den
    return [sign,n]


def model():
    ctx.prec,ctx.cap=480,3
    _,_,gs,_,weights=joint["joint_data"]()
    nums=[tick(w) for w in weights[:-1]]
    nums.append(den-sum(nums))
    rows=[]
    for i,w,n in zip(joint["support"],weights,nums):
        g=gs[i-1]
        assert abs(arb(n)/den-w)<arb(8)/den
        rows.append([i,tick(g),n,rational_turn(2*g*arb(2).log()),rational_turn(2*g*arb(3).log())])
    return {"bits":bits,"rows":rows,"gram_bins":original["bins"]}


def main():
    data=json.loads((directory/"ch39-integer-shell-record.json").read_text())
    assert data==model()
    print("Integer heights, both rational turns and the weights calibrated at all eight actual ordinates.")


if __name__=="__main__":
    main()

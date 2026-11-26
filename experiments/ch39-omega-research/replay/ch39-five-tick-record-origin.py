"""Calibrate the five-tick rational record using the 46 certified zeta heights."""
import json
from pathlib import Path
from runpy import run_path

from flint import arb, acb, arb_mat, ctx

directory=Path(__file__).parent
setup=run_path(str(directory/"ch39-rh-shell-smallest-tick.py"))
single=run_path(str(directory/"ch39-integer-shell-record-origin.py"))
bits=288
den=2**bits


def code(v):
    k=(v*den).floor().unique_fmpz()
    assert k is not None
    return int(k)


def model():
    ctx.prec,ctx.cap=1664,3
    levels=setup["prior"]["read"](42)
    heights=[acb.zeta_zero(i).imag for i in setup["support"]]
    logs=[5*arb(t).log() for t in (2,3,arb(3)/2,6)]
    M=arb_mat([[g*g/(j*j+g*g) for j in range(42)]+[(g*t).cos() for t in logs] for g in heights])
    weights=M.transpose().solve(arb_mat([[v] for v in levels]+[[0]]*4)).entries()
    assert min(weights)>arb(".00079")
    return levels,heights,weights


def chart(theta):
    x,y=theta.cos(),theta.sin()
    s=1 if x>0 else -1
    k=code(s*y/(1+s*x))
    return [s,k]


def record():
    levels,heights,weights=model()
    nums=[code(w) for w in weights[:-1]]
    nums.append(den-sum(nums))
    rows=[[i,code(g),n,chart(5*g*arb(2).log()),chart(5*g*arb(3).log())]
          for i,g,n in zip(setup["support"],heights,nums)]
    targets=[code(f) for f in levels]
    return {"bits":bits,"levels":targets,"rows":rows,"gram_bins":single["original"]["bins"]}


if __name__=="__main__":
    assert record()==json.loads((directory/"ch39-five-tick-record.json").read_text())
    print("42 dyadic RH levels, 46 actual-height bins, integer counts and two rational turns per row certified.")

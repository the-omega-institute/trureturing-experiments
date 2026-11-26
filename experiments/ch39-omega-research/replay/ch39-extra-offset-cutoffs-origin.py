"""Calibrate or verify the rational support and separation witnesses at unit scale."""
from pathlib import Path
from runpy import run_path
import json

from flint import acb, arb, arb_mat, arb_series, ctx


def level(j, mass):
    if j==0:
        return arb(1)
    t=arb_series([arb(1)/2+j,1])
    xi=t*(t-1)/2*((-t/2)*arb.pi().log()).exp()*(t/2).gamma()*t.zeta()
    return xi[1]/(2*j*xi[0]*mass)


def compute(case):
    root=Path(__file__).resolve().parent
    words=run_path(str(root/"ch39-shortest-triangle-clock-words.py"))
    rates=run_path(str(root/"ch39-two-prime-height-law.py"))["rate_numerators"]()
    ctx.prec,ctx.cap=1800,3
    t=arb_series([arb(1)/2,1])
    xi=t*(t-1)/2*((-t/2)*arb.pi().log()).exp()*(t/2).gamma()*t.zeta()
    mass=xi[2]/xi[0]
    offsets=case["offsets"]; ids=case["support"]; size=len(ids)
    gs=[acb.zeta_zero(i).imag for i in ids]
    angular=[[] for _ in range(4)]
    for g in gs:
        h=int((g*2**288).floor().unique_fmpz())
        m,n=[-43235114615559591407865398948523*((h*r+2**311)//2**312) for r in rates]
        for row,k in zip(angular,(m,n,n-m,n+m)):
            lo,hi=words["power"](k)[0]
            row.append((arb(lo)/words["DEN"]).union(arb(hi)/words["DEN"]))
    rows=[[1/(1+(j/g)**2) for g in gs] for j in offsets]+angular
    assert size==len(rows) and max(ids)==case["cut"]
    inv=arb_mat([[v.mid() for v in row] for row in rows]).inv()
    ident=arb_mat([[int(i==j) for j in range(size)] for i in range(size)])
    E=inv*arb_mat(rows)-ident
    rho=max(sum(abs(E[i,j]).upper() for j in range(size)).upper() for i in range(size))
    assert rho<arb(2)**-280
    trial=(inv*arb_mat([[level(j,mass)] for j in offsets]+[[0]]*4)).entries()
    pad=rho*max(abs(v).upper() for v in trial)/(1-rho)
    assert min(v.lower() for v in trial)-pad>arb(".00004")
    counts=[int((v.mid()*2**288).floor().unique_fmpz()) for v in trial[:-1]]
    counts.append(2**288-sum(counts))
    anchors=[1/acb.zeta_zero(i).imag**2 if i else arb(0) for i in case["anchors"]]
    M=arb_mat([[1/(1+j*j*z) for j in offsets] for z in anchors])
    coef=M.solve(arb_mat([[arb(v)/2**48] for v in case["dual"][:-4]])).entries()
    a=[int((v*2**40).floor().unique_fmpz()) for v in coef]
    return a,counts


if __name__=="__main__":
    root=Path(__file__).resolve().parent
    cases=json.loads((root/"ch39-extra-offset-cutoffs.json").read_text())
    for c in cases:
        a,weights=compute(c)
        assert c["a"]==a and c["weights"]==weights
    print("Both exact positive realizations and rational separating rows verified.")


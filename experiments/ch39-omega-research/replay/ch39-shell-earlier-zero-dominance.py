"""Certified earlier-zero dominance for geometric energies monotone in the four clock strengths."""
from fractions import Fraction as F
from itertools import product
from flint import arb,acb,ctx
ctx.prec=96
freq=[]
for r,d0,d1,d2 in ((17,6,12,14),(867,244,724,804)):
    theta=[(1-arb(d)/(2*r)).acos() for d in (d0,d1,d2)]
    freq.extend([x-theta[0] for x in theta[1:]])
cs=[(arb(x.numerator)/x.denominator).sqrt() for x in (F(81,98),F(8741547,12563572))]
def obs(h):
    phases=[(h*f).sin_cos() for f in freq]
    strength=[4*s*s*c*c for s,c in phases]
    triangle=[]
    for i,c in zip((0,2),cs):
        (s,u),(t,v)=phases[i:i+2]
        triangle.extend((s*s,t*t,1-(u*v+c*s*t)**2))
    return strength,triangle
target=obs(acb.zeta_zero(719).imag)
cells=[{},{}]
both={"lower":[],"higher":[]}
for start in range(1,719,64):
    for n,z in enumerate(acb.zeta_zeros(start,min(64,719-start)),start=start):
        actual=obs(z.imag)
        signs=[]
        for row,ref,table in zip(actual,target,cells):
            diffs=[u-v for u,v in zip(row,ref)]
            assert all(not z.contains(0) for z in diffs)
            mask=tuple(-1 if z<0 else 1 for z in diffs)
            table.setdefault(mask,[]).append(n)
            signs+=list(mask)
        if max(signs)<0:
            both["lower"].append(n)
        if min(signs)>0:
            both["higher"].append(n)
for table in cells:
    width=len(next(iter(table)))
    print("observed orthants of",width,"dimensions",len(table),"of",2**width)
    for mask in [(-1,)*width,(1,)*width]:
        ids=table.get(mask,[])
        print("earlier zeros with all signs",mask[0],len(ids),"first",ids[:5])
    if len(table)<2**width:
        print("missing sign patterns",sorted(set(product((-1,1),repeat=width))-set(table)))
print("simultaneously dominating all ten quantities",both)

assert len(cells[0])==16 and len(cells[1])==46 and both["lower"]==[62]

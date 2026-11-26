"""Spin-one clock axes as family adjoints and a dimension-seven Majorana test texture."""
from pathlib import Path
from fractions import Fraction as F
from flint import acb,arb,fmpz_mat
import runpy

book=runpy.run_path(str(Path(__file__).with_name("ch39-book-shell-clock-embedding.py")))
def flatten(x):
    return [x[i,j] for i in range(3) for j in range(3)]
I=fmpz_mat([[1,0,0],[0,1,0],[0,0,1]])
# A common factor i does not affect the complex span of the words at any length.
X=fmpz_mat([[0,0,0],[0,0,-1],[0,1,0]])
Z=fmpz_mat([[0,-1,0],[1,0,0],[0,0,0]])
words=[I]
all_words=words[:]
ranks=[1]
for depth in range(1,4):
    words=[w*g for w in words for g in (X,Z)]
    all_words+=words
    ranks.append(fmpz_mat([flatten(x) for x in all_words]).rank())
assert ranks==[1,3,7,9]
for r,offset,c2 in ((17,0,F(81,98)),(867,2,F(8741547,12563572))):
    c2=arb(c2.numerator)/c2.denominator
    p,q=[z*book["h"] for z in book["f"][offset:offset+2]]
    a,b=(2*p).sin(),(2*q).sin()
    assert not (a*b).contains(0)
    trace=a*a+b*b
    disc=((a*a-b*b)**2+4*a*a*b*b*c2).sqrt()
    high=(trace+disc)/2
    low= a*a*b*b*(1-c2)/high
    assert low>0 and high>low and (low>1 or low<1) and (high>1 or high<1)
    masses=sorted([arb(1),1/high,1/low],key=lambda v:float(v.mid()))
    gaps=[masses[i+1]**2-masses[i]**2 for i in (0,1)]
    assert min(gaps)>0
    ratio=min(gaps)/sum(gaps)
    print("squared shell radius",r,"heavy Takagi ratios 1,",high.str(13),low.str(13))
    print("seesaw light masses in common units",*[z.str(12) for z in masses],
          "small adjacent squared gap / total spread",ratio.str(13))

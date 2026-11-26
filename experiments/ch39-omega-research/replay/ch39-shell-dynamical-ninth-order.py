"""Release the ninth-order reference angles and minimize over them together with the zero."""
from pathlib import Path
import runpy
from flint import arb,acb

data=runpy.run_path(str(Path(__file__).with_name("ch39-shell-earlier-zero-dominance.py")))
def triples(c):
    r=[]
    for i in range(9):
        for j in range(9):
            if (i,j)>((-i)%9,(-j)%9):
                continue
            a,b=(arb.pi()*i/9).sin_cos()
            t,v=(arb.pi()*j/9).sin_cos()
            r.append(((i,j),(a*a,t*t,1-(b*v+c*a*t)**2)))
    return r
tables=[triples(c) for c in data["cs"]]
lam=arb(1)/51**5
def distance(x,y):
    return sum((a-b)**2 for a,b in zip(x,y))
def best(x,table):
    rows=sorted(((distance(x,r),t,r) for t,r in table),key=lambda v:float(v[0].mid()))
    # Some isolated heat angles may distinguish ties by less than numeric precision;
    # this computation is on a short finite list and verifies the winning interval.
    assert all(rows[0][0]<v[0] for v in rows[1:])
    return rows[0]
star=[dict(t)[v] for t,v in zip(tables,((2,1),(4,6)))]
T=acb.zeta_zero(719).imag
old=data["obs"](T)[1]
estar=distance(old,star[0]+star[1])+lam*T*T
champion=None
scores={}
dominance=[]
n=1
while True:
    h=acb.zeta_zero(n).imag
    if champion and lam*h*h>champion[0]:
        break
    x=data["obs"](h)[1]
    a,b=[best(xx,t) for xx,t in zip((x[:3],x[3:]),tables)]
    score=a[0]+b[0]+lam*h*h
    scores[n]=score
    if champion is None or score<champion[0]:
        champion=(score,n,a[1]+b[1])
    elif n!=champion[1]:
        assert score>champion[0]
    if score<estar and all(u<=v or u.overlaps(v) for u,v in zip(a[2]+b[2],star[0]+star[1])):
        dominance.append((n,a[1]+b[1],score))
    n+=1
assert all(champion[0]<v for k,v in scores.items() if k!=champion[1])
assert champion[1]==98 and champion[2]==(4,0,4,5)
print("unfixed reference ground state",champion[1],champion[2],champion[0].str(15),
      "first height excluded by floor",n)
print("first weakly dominating references",[(j,s,v.str(10)) for j,s,v in dominance[:3]])
# Search for a strict six-coordinate dominance witness at an earlier height; only this
# strict test is used for the statement on arbitrary monotone reference penalties.
for j in range(1,720):
    h=acb.zeta_zero(j).imag
    if lam*h*h>=estar:
        break
    x=data["obs"](h)[1]
    local=[]
    for block,xx,ss in zip(tables,(x[:3],x[3:]),star):
        vals=[(distance(xx,r),idx,r) for idx,r in block if all(u<v for u,v in zip(r,ss))]
        assert vals
        local.append(min(vals,key=lambda z:float(z[0].mid())))
    total=sum(z[0] for z in local)+lam*h*h
    if total<estar:
        assert j==62
        print("strict dominance witness",j,local[0][1]+local[1][1],total.str(15),"719 score",estar.str(15))
        break
else:
    raise AssertionError("strict dominance witness not found")

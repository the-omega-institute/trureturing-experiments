"""Four positive pentagonal actions expose the fitted clock pair on the full torus."""
from pathlib import Path
import runpy
from flint import arb,arb_mat

loop = runpy.run_path(str(Path(__file__).with_name("ch39-clock-wilson-length.py")))
m,q = loop["m"],loop["q"]
pi = arb.pi()
axis = m["m"]["c"]
paths = [(0,6,8,3,7),(0,6,1,9,8),(0,6,3,7,8),(0,6,7,2,8)]
powers = [(2,0),(1,1),(0,2),(1,-2)]

def phases():
    reals = []
    for path,(k,l) in zip(paths,powers):
        phase=q(1)
        word=[]
        for edge in zip(path,path[1:]+path[:1]):
            pair=tuple(sorted(edge))
            d,z = loop["links"][pair]
            phase*=z if pair==edge else z.conjugate()
            if d in (724,804):
                word.append((1 if d==724 else 2)*(1 if pair==edge else -1))
        word=loop["key"](word)
        # Both letters can be gathered into one block of each type by cyclicity.
        assert len(word)<=3 and (sum(abs(i)==1 for i in word),sum(abs(i)==2 for i in word))==(abs(k),abs(l))
        pair = (word.count(1)-word.count(-1),word.count(2)-word.count(-2))
        assert pair in ((k,l),(-k,-l))
        reals.append(phase.real)
    return reals

berry = phases()
_,angles = m["rotations"](719)
target=[x-2*pi*round(float((x/(2*pi)).mid())) for x in angles]
assert target[0]>0

def data(p,t,order):
    rows = []
    for (k,l),z in zip(powers,berry):
        # t = (1+c) cos(kp+lq)+(1-c) cos(kp-lq).
        x=k*p+l*t
        y=k*p-l*t
        if order==0:
            u=(1+axis)*x.cos()
            v=(1-axis)*y.cos()
        else:
            sx,cx=x.sin_cos()
            sy,cy=y.sin_cos()
            u,v=(1+axis)*cx,(1-axis)*cy
            a,b=(1+axis)*sx,(1-axis)*sy
        tr=u+v
        result=[8-tr*tr-2*z*tr]
        if order:
            tp,tt=-k*(a+b),-l*(a-b)
            result.extend((-2*(tr+z)*tp,-2*(tr+z)*tt))
        if order==2:
            result.extend((2*((tr+z)*k*k*tr-tp*tp),
                           2*((tr+z)*k*l*(u-v)-tp*tt),
                           2*((tr+z)*l*l*tr-tt*tt)))
        rows.append(result)
    return rows

station=data(*target,1)
C=arb_mat([[station[j][i] for j in (1,2)] for i in (1,2)])
fixed=q("1.8")
co=C.solve(arb_mat([[-station[0][i]-fixed*station[3][i]] for i in (1,2)]))
weights=[q(1),co[0,0],co[1,0],fixed]
assert min(w.lower() for w in weights)>0
def eval_cost(p,t,order):
    return [sum(c*z for c,z in zip(weights,col)) for col in zip(*data(p,t,order))]

emin=eval_cost(*target,0)[0]
radius=q(".10")
def local_check():
    lower=q(5)
    for i in range(12):
        p=target[0]+radius*(q(i).union(q(i+1))/6-1)
        for j in range(12):
            t=target[1]+radius*(q(j).union(q(j+1))/6-1)
            a,b,c=eval_cost(p,t,2)[3:]
            assert a>lower and (a-lower)*(c-lower)-b*b>0
    return lower

def exterior():
    # Simultaneous inversion is exact, so the half-torus covers both wells.
    work=[(q(0).union(pi),(-pi).union(pi),0)]
    count=0
    while work:
        p,t,d=work.pop()
        if abs(p-target[0])<radius and abs(t-target[1])<radius:
            continue
        if eval_cost(p,t,0)[0]-emin>q(".015"):
            count+=1
            continue
        assert d<32
        x = p if d%2==0 else t
        mid=x.mid()
        for y in (x.lower().union(mid),mid.union(x.upper())):
            work.append((y,t,d+1) if d%2==0 else (p,y,d+1))
    return count

if __name__=="__main__":
    local=local_check()
    count=exterior()
    print("weights",*[w.str(11) for w in weights],"local Hessian lower",local)
    print("positive exterior cells",count,"stationarity matrix determinant",C.det().str(12))

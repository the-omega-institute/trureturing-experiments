"""所有实现同一十四项权的字概率中，最优公平停时上限为 3924。"""

from itertools import product
from math import comb
from pathlib import Path
from runpy import run_path
import sys

import numpy as np
from flint import arb, arb_mat, ctx


WINNER = (1,1,1,-1,1,-1,-1,-1,-1,-1,-1,1,1,1)
OTHER = (1,1,1,-1,1,-1,-1,-1,-1,1,-1,1,1,1)
MULT = [1]+[2]*106


def planes(root):
    load = lambda name: run_path(str(root/("ch39-"+name+".py")))
    shell = load("three-ray-chiral-boundary")
    model = load("mincut-global-winding")
    local = load("mincut-reference-ckm-action")
    flavor = load("chiral-neutral-flavor-phase")
    ctx.prec = 256
    x = arb("-.01046")
    y = x-local["reference_angles"](flavor)[0]
    weights = load("height-independent-flavor-vacuum")["minimum_order_weights"](
        model,model["terms"](shell),x,y)
    support = load("realized-vacuum-word-masses")["SUPPORT"]
    theta = load("mincut-real-triangle-tomography")["unit"]((108469,109)).arg()
    matrix = arb_mat([[(q*j*theta).cos() for j in range(107)] for q in range(107)])
    inverse = matrix.inv()
    base = [inverse[j,0] for j in range(107)]
    free = [inverse[j,106] for j in range(107)]
    columns = [[inverse[j,q]*weights[q].sqrt() for j in range(107)] for q in support]
    return matrix, weights, base, free, columns


def pairs(base, free, columns):
    # 浮点只建议两个非负约束作为对偶支点。
    a,b = [np.array(list(map(float,l))) for l in (base,free)]
    basis = np.array([[float(x) for x in row] for row in columns])
    j = np.arange(107)
    ca,cb = a@j,b@j
    i,k = np.triu_indices(107,1)
    signs = np.array([(1,*tail) for tail in product((1,-1),repeat=13)])
    for batch in np.array_split(signs,128):
        v = batch@basis
        cv = (v@j)[:,None]
        determinant = v[:,i]*b[k]-v[:,k]*b[i]
        with np.errstate(divide="ignore",invalid="ignore"):
            alpha = (cv*b[k]-cb*v[:,k])/determinant
            beta = (cb*v[:,i]-cv*b[i])/determinant
        values = np.where((alpha>0)&(beta>0),ca-alpha*a[i]-beta*a[k],-np.inf)
        for sign,index in zip(batch,values.argmax(axis=1)):
            yield tuple(map(int,sign)), int(i[index]), int(k[index])


def abs_cost(vector):
    return sum(j*x for j,x in enumerate(vector))


def lower_bound(a, v, b, cost, i, j):
    ca,cv,cb = cost
    det = v[i]*b[j]-v[j]*b[i]
    alpha = (cv*b[j]-cb*v[j])/det
    beta = (cb*v[i]-cv*b[i])/det
    assert not det.contains(0) and alpha>0 and beta>0
    return ca-alpha*a[i]-beta*a[j]


def calls(groups):
    return [sum(((j-k)*groups[j] for j in range(k+1,107)),arb(0)) for k in range(107)]


def affine_obstruction(a,v,b,probe):
    # 浮点选择每一步在哪些位置保留仿射分支。选择哪一支都给严格有效的下界。
    p = np.array(list(map(float,probe)))
    j = np.arange(107)
    guide = np.maximum(j[None,:]-j[:,None],0)@p
    probability = p/np.array(MULT)
    affine = [calls(row) for row in (a,v,b)]
    masses = [[x/n for x,n in zip(row,MULT)] for row in (a,v,b)]
    zero = arb(0)
    for _ in range(3923):
        guide = np.maximum(0,np.r_[guide[1],(guide[:-2]+guide[2:])/2,0]-probability)
        guide[-1] = 0
        mask = np.flatnonzero(guide>0).tolist()
        for k,(row,mass) in enumerate(zip(affine,masses)):
            following = [zero]*107
            for x in mask:
                following[x] = (row[abs(x-1)]+row[x+1])/2-mass[x]
            affine[k] = following
    return lower_bound(a,v,b,[row[1] for row in affine],1,50)


def upper(groups):
    zero = arb(0)
    mu = [x/d for x,d in zip(groups,MULT)]
    assert min(mu)>arb("4.7e-7") and (sum(groups)-1).contains(0)
    gap = calls(groups)
    active = [arb(1)]+[zero]*106
    held = [zero]*107
    for step in range(3924):
        following = [zero]*107
        released = [zero]*107
        for j,d in enumerate(gap):
            if d.is_zero():
                held[j]+=active[j]
                continue
            assert j<106 and held[j].is_zero()
            candidate = (gap[abs(j-1)]+gap[j+1])/2-mu[j]
            if candidate>0:
                following[j]=candidate
                released[j]=active[j]
            else:
                assert candidate<0 and active[j]>d
                released[j]=d
                held[j]+=active[j]-d
        if step==3923:
            assert all(z.is_zero() for j,z in enumerate(gap) if j!=1)
            assert gap[1]>0
            print("正概率见证最后一步的总放走量",(2*gap[1]).str(17))
        assert released[106].is_zero()
        active = [released[1]]+[(x+y)/2 for x,y in zip(released[:-2],released[2:])]+[released[105]/2]
        gap = following
    assert all(z.is_zero() for z in gap)
    errors = [h+w-p for h,w,p in zip(held,active,mu)]
    assert all(e.contains(0) and abs(e)<arb("1e-35") for e in errors)
    print("见证平均步数、每个字的最小概率",
          sum(k*k*p for k,p in enumerate(groups)).str(18),min(mu).str(16))


def verify(root):
    matrix,weights,a,b,columns = planes(root)
    threshold = arb(3923*comb(3922,1961))/2**3922
    ca,cb = abs_cost(a),abs_cost(b)
    cs = [abs_cost(row) for row in columns]
    skipped = 0
    other_floor = arb(100)
    minimum = arb(100)
    for signs,i,j in pairs(a,b,columns):
        cv = sum(h*z for h,z in zip(signs,cs))
        v = {k:sum(h*row[k] for h,row in zip(signs,columns)) for k in (i,j)}
        bound = lower_bound(a,v,b,(ca,cv,cb),i,j)
        minimum = minimum.min(bound)
        if signs not in (WINNER,OTHER):
            assert bound>threshold
            skipped+=1
            other_floor = other_floor.min(bound)
            continue
        assert (i,j)==(1,50)
        v = [sum(h*row[k] for h,row in zip(signs,columns)) for k in range(107)]
        det=v[i]*b[j]-v[j]*b[i]
        s=(a[j]*b[i]-a[i]*b[j])/det
        t=(a[i]*v[j]-a[j]*v[i])/det
        probe=[x+s*y+t*z for x,y,z in zip(a,v,b)]
        assert probe[i].contains(0) and probe[j].contains(0)
        probe[i]=probe[j]=arb(0)
        assert min(p for k,p in enumerate(probe) if k not in (i,j))>0
        assert (abs_cost(probe)-bound).contains(0)
        residual = affine_obstruction(a,v,b,probe)
        assert residual>arb(".00110" if signs==WINNER else ".04751")
        print("剩余平面的符号、最小绝对位置均值、3923时位置1的统一势差下界",
              signs,bound.str(17),residual.str(17))
        if signs==WINNER:
            eps=arb(1)/10000
            witness=[(1-eps)*p+eps*x for p,x in zip(probe,a)]
            moments=(matrix*arb_mat([[p] for p in witness])).entries()
            norm=sum(z*z for z in moments[1:106])
            assert norm>0 and all((moments[q]**2/norm-weights.get(q,0)).contains(0) for q in range(1,106))
            assert (moments[0]-1).contains(0)
            upper(witness)
    assert skipped==8190 and minimum>arb("49.90098186280") and other_floor>arb("49.99323")
    print("全局最小绝对位置均值、3923 步公平方向的绝对位置均值、已排除平面数及共同下界",
          minimum.str(17),threshold.str(17),skipped,other_floor.str(17))


if __name__=="__main__":
    root=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
    verify(root)

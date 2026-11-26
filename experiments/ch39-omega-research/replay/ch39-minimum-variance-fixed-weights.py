"""固定十四项绕数权，认证所有共轭对称字概率中最小的二阶矩。"""

from itertools import product
from pathlib import Path
from runpy import run_path
import sys

import numpy as np
from flint import arb, arb_mat, ctx


OPTIMAL_SIGNS = (1,-1,1,-1,-1,1,-1,-1,1,1,1,1,1,1)


def second_moment(row):
    return sum(j*j*p for j,p in enumerate(row))


def proposed_pairs(base, free, columns):
    # 给出二维线性规划的两个对偶支点。所有符号、乘子及下界另以 Arb 验证。
    a, b = (np.array(list(map(float, row))) for row in (base, free))
    matrix = np.array([[float(x) for x in row] for row in columns])
    objective = np.arange(107)**2
    cost_a, cost_b = a@objective, b@objective
    i,j = np.triu_indices(107,1)
    signs = np.array([(1,*tail) for tail in product((1,-1),repeat=len(columns)-1)])
    for batch in np.array_split(signs,128):
        v = batch@matrix
        cost_v = (v@objective)[:,None]
        det = v[:,i]*b[j]-v[:,j]*b[i]
        with np.errstate(divide="ignore", invalid="ignore"):
            left = (cost_v*b[j]-cost_b*v[:,j])/det
            right = (cost_b*v[:,i]-cost_v*b[i])/det
        score = np.where((left>0)&(right>0),cost_a-left*a[i]-right*a[j],-np.inf)
        for sign,k in zip(batch,np.argmax(score,axis=1)):
            yield tuple(map(int,sign)), int(i[k]), int(j[k])


def verify(root):
    load = lambda name: run_path(str(root/("ch39-"+name+".py")))
    published = load("realized-vacuum-word-masses")
    target = load("height-independent-flavor-vacuum")
    model = load("mincut-global-winding")
    geometry = load("mincut-real-triangle-tomography")
    local = load("mincut-reference-ckm-action")
    shell = load("three-ray-chiral-boundary")
    flavor = load("chiral-neutral-flavor-phase")
    ctx.prec = 224
    x = arb("-.01046")
    y = x-local["reference_angles"](flavor)[0]
    weights = target["minimum_order_weights"](model,model["terms"](shell),x,y)
    support = published["SUPPORT"]
    assert tuple(sorted(weights)) == support and len(support)==14
    theta = geometry["unit"]((108469,109)).arg()
    matrix = arb_mat([[(q*j*theta).cos() for j in range(107)] for q in range(107)])
    inverse = matrix.inv()
    base = [inverse[j,0] for j in range(107)]
    free = [inverse[j,106] for j in range(107)]
    columns = [[inverse[j,q]*weights[q].sqrt() for j in range(107)] for q in support]
    assert min(base)>0
    cost_a,cost_b = map(second_moment,(base,free))
    cost_columns = [second_moment(c) for c in columns]
    runner_up = arb(10**6)
    planes = 0
    for signs,i,j in proposed_pairs(base,free,columns):
        planes += 1
        vi,vj = [sum(h*c[k] for h,c in zip(signs,columns)) for k in (i,j)]
        cost_v = sum(h*z for h,z in zip(signs,cost_columns))
        det = vi*free[j]-vj*free[i]
        assert not det.contains(0)
        left = (cost_v*free[j]-cost_b*vj)/det
        right = (cost_b*vi-cost_v*free[i])/det
        assert left>0 and right>0
        lower = cost_a-left*base[i]-right*base[j]
        if signs != OPTIMAL_SIGNS:
            runner_up = runner_up.min(lower)
            continue
        assert (i,j)==(17,60)
        scale = (base[j]*free[i]-base[i]*free[j])/det
        tail = (base[i]*vj-base[j]*vi)/det
        direction = [sum(h*c[k] for h,c in zip(signs,columns)) for k in range(107)]
        groups = [a+scale*v+tail*b for a,v,b in zip(base,direction,free)]
        assert groups[i].contains(0) and groups[j].contains(0)
        # 这两项由定义 s,t 的非奇异线性方程精确为零。
        groups[i] = groups[j] = arb(0)
        assert min(g for k,g in enumerate(groups) if k not in (i,j))>arb(".000238")
        assert scale>0 and (sum(groups)-1).contains(0)
        value = second_moment(groups)
        assert (value-lower).contains(0)
    assert planes==8192 and runner_up>arb("3381.23319")
    assert arb("3379.59587479924") < value < arb("3379.59587479926")
    moments = (matrix*arb_mat([[p] for p in groups])).entries()
    assert (moments[0]-1).contains(0) and (moments[106]-tail).contains(0)
    total = sum(z*z for z in moments[1:106])
    assert total>0 and (total-scale**2).contains(0)
    assert all((moments[q]**2/total-weights.get(q,0)).contains(0) for q in range(1,106))
    previous = load("fair-vacuum-horizon")["target_masses"](root)
    saving = second_moment(previous)-value
    tv = sum(abs(z-arb(n)/213) for z,n in zip(groups,[1]+[2]*106))/2
    print("符号平面数、其余平面二阶矩共同下界",planes,runner_up.lower().str(18))
    print("最小二阶矩、减少的平均步数、距逐字均匀的总变差",
          *(z.str(19) for z in (value,saving,tv)))
    print("矩的尺度、第106矩、取零的 |j|",scale.str(19),tail.str(19),17,60)
    print("按 |j|=0,…,106 列每个单字的最优概率")
    for start in range(0,107,10):
        print(start,*[(groups[k]/(1 if k==0 else 2)).str(13) for k in range(start,min(start+10,107))])


if __name__=="__main__":
    root = Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
    verify(root)

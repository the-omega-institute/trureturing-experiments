"""两背景各自达到两引力矩时，连续共同插值分布在大标尺上的规范障碍。"""
from concurrent.futures import ProcessPoolExecutor
from functools import cache
import json
from os import cpu_count
from math import factorial
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

ROOT = Path(__file__).resolve().parent

# 来自有限格点线性规划的有理候选；下面在两维连续域独立验证其余量。
DIRECTIONS = json.loads(Path(__file__).with_name("ch39-shared-large-scale-directions.json").read_text())


def load():
    earlier = run_path(str(ROOT/'ch39-shared-continuous-position-gap.py'))
    source = earlier['load']('shared-twelve-node-gauge-gap')
    kernel = earlier['load']('continuous-split-roots')
    _,lib,backgrounds = earlier['load']('two-zero-cutoff-gravity')['gravity_fields'](
        ROOT,source['HEAVY'],tight_neutral=True)
    model,_,_,initial,_=lib['chiral']['source']['mean']['inputs']()
    fields=[(0,0,0)]+earlier['load']('signed-gravity-unification')['gauge_weights'](model)
    shell=earlier['load']('native-rh-ray-pvm')['shell_weights']
    beta=[v/(2*arb.pi()) for v in model['beta']]
    limit=(arb('1.22089e19')/model['mz']).log()
    data=[]
    for n,bg in enumerate(backgrounds):
        edges=source['paired_lines'](ROOT,shell,source['SHELL'][n])
        assert len(edges)==361 and all(0 < a and (a<b or a==b) and b<=1 for a,b,_ in edges)
        d=9*min(a.lower() for a,b,w in edges)/1000
        assert max(z.upper() for g,k,z in bg)<min(arb(1),5*d)
        data.append((edges,kernel['series'](bg,fields),bg))
    return data,fields,beta,limit,initial,kernel,earlier


def project_direction(ctx, components, kernel_ratio=5, *, low_colours=None, common_log=None):
    data,fields,beta,limit,initial,kernel,earlier=ctx
    y=[arb(v)/10**11 for v in components]
    assert len(y)==4*len(data)
    norm=sum(abs(v) for j,v in enumerate(y) if j%4 in (2,3))
    tables=[]
    constant=arb(0)
    order=kernel['ORDER']
    assert arb(kernel_ratio)/(order+2) < arb(1)/2
    tail=2*arb(kernel_ratio)**(order+1)/factorial(order+1)
    for n,(edges,series,bg) in enumerate(data):
        ky,wy,u,v=y[4*n:4*n+4]
        colour = arb(0) if low_colours is None else low_colours[n]
        direction=[ky,wy,u,v,colour-u-v]
        row=arb_mat([direction])*series
        tables.append([(a,b-a,w*row) for a,b,w in edges])
        rate=u*(beta[0]-beta[2])+v*(beta[1]-beta[2])
        shift = (rate*limit).max(0) if common_log is None else rate*common_log
        constant+= -ky-wy*arb('1e-122')+u*(initial[0]-initial[2])+v*(initial[1]-initial[2])+colour*initial[2]-shift
        size=sum(abs(2*ky*k)+abs(wy*g)+sum(abs(d*j) for d,j in zip(direction[2:],bs,strict=True))
                 for (g,k,z),bs in zip(bg,fields,strict=True))/(4*arb.pi())
        constant-=tail*size*sum(abs(w) for a,b,w in edges)
    assert norm > 0
    return list(zip(*tables,strict=True)),constant,norm


def curvature(rows,unit,lo):
    """c 的二阶导数。两背景二次主项逐行合并；对数项跨行求和。"""
    logarithm=arb(0)
    upper=arb(0)
    for pair in rows:
        b0,b1,b2=arb(0),arb(0),arb(0)
        for a,span,coeff in pair:
            b0+=coeff[0,0]*(a/unit)**2
            b1+=coeff[0,0]*(a/unit)*((a+span)/unit)
            b2+=coeff[0,0]*((a+span)/unit)**2
            logarithm+=coeff[0,coeff.ncols()-1]
            ratio=unit/(lo*a)
            power=arb(1)
            for j in range(1,coeff.ncols()-3):
                power*=ratio
                upper+=j*(j+1)*abs(coeff[0,j+2])*power/lo**2
        # 二次 Bernstein 系数控制任意位置分布的平均二阶导数上界。
        upper+=2*max(b.upper() for b in (b0,b1,b2))
    # 下界只需 f'' 的上界。负的同向项也计入这个绝对上界。
    return (upper+abs(logarithm)/lo**2).max(0).upper()


def at(ctx,model,micro):
    kernel=ctx[-2]
    rows,constant,norm=model
    c=arb(micro)/10**6
    minimum=ctx[-1]['projected_minimum']
    lower=constant
    for pair in rows:
        lower+=minimum(kernel['basis'],kernel['S'],[(c*a,c*span,row) for a,span,row in pair],arb('0.0001'))
    return lower.lower()


def certify_part(ctx,key,left,right):
    model=project_direction(ctx,DIRECTIONS[str(key)])
    rows,_,norm=model
    M=curvature(rows,ctx[-2]['S'],arb(left)/10**6)
    @cache
    def endpoint(n):
        return at(ctx,model,n)
    todo=[(left,right)]
    completed=[]
    while todo:
        a,b=todo.pop()
        l,r=endpoint(a),endpoint(b)
        width=arb(b-a)/10**6
        bound=(min(l,r)-M*width**2/8)/norm
        if bound > arb('.05'):
            completed.append(bound.lower())
        else:
            assert min(l,r)>arb('.05')*norm, (key,a,b,l,r,norm)
            m=(a+b)//2
            assert a<m<b
            todo.extend(((a,m),(m,b)))
    return min(completed),len(completed)


@cache
def inputs():
    return load()


def certify_segment(interval):
    key,left,right=interval
    bound,cells=certify_part(inputs(),key,left,right)
    return bound.str(20),cells


def verify():
    keys=sorted(map(int,DIRECTIONS))
    assert keys[0]==9000 and keys[-1]==1000000
    cuts=[keys[0],*[(a+b)//2 for a,b in zip(keys,keys[1:])],keys[-1]]
    parts=list(zip(keys,cuts[:-1],cuts[1:],strict=True))
    with ProcessPoolExecutor(max_workers=min(8,cpu_count() or 1)) as pool:
        results=list(pool.map(certify_segment,parts))
    bound=min(arb(v).lower() for v,n in results)
    assert bound > arb('.06354')
    print('标尺区间 [0.009,1]、分离方向数、完整闭段数、规范差下界',
          len(parts),sum(n for _,n in results),bound.str(12),flush=True)


if __name__=='__main__':
    verify()

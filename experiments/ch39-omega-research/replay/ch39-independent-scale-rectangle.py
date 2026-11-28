"""在独立标尺二维闭矩形上，认证任意连续逐线共同分布的联合交汇障碍。"""
from concurrent.futures import ProcessPoolExecutor
from functools import cache, lru_cache
from itertools import product
import json
from os import cpu_count
from pathlib import Path
from runpy import run_path

from flint import arb

ROOT = Path(__file__).resolve().parent
UNIT = 10**6
TARGET = arb('0.01')
BOX = (9500, 10500, 7500, 9000)
DIRECTIONS = json.loads((ROOT/'ch39-independent-scale-directions.json').read_text())


@cache
def inputs():
    high = run_path(str(ROOT/'ch39-shared-large-scale-gap.py'))
    ctx = high['load']()
    for ix,(edges, _, bg) in enumerate(ctx[0]):
        floor = arb(BOX[2*ix])/UNIT*min(a.lower() for a,_,_ in edges)
        assert max(z.upper() for _,_,z in bg) < min(arb(1), 6*floor)
    return high, ctx


@lru_cache(maxsize=8)
def prepare(key):
    high, ctx = inputs()
    rows, const, norm = high['project_direction'](ctx,DIRECTIONS[key],6)
    curvature = [high['curvature']([(pair[ix],) for pair in rows],
                                    ctx[-2]['S'],arb(BOX[2*ix])/UNIT) for ix in (0,1)]
    return rows, const, norm, curvature


@cache
def endpoint(key, micro0, micro1):
    _, ctx = inputs()
    rows,lower,norm,_ = prepare(key)
    scales = [arb(x)/UNIT for x in (micro0,micro1)]
    for pair in rows:
        line = [(c*a,c*d,coeff) for c,(a,d,coeff) in zip(scales,pair,strict=True)]
        lower += ctx[-1]['projected_minimum'](ctx[-2]['basis'],ctx[-2]['S'],line,arb('1e-9'))
    return (lower/norm).lower()


def attempt(key, box):
    x0,x1,y0,y1=box
    _,_,norm,(mx,my)=prepare(key)
    cost=(mx*(x1-x0)**2+my*(y1-y0)**2)/(8*UNIT**2*norm)
    # 每个背景只用各自标尺，角点的双线性插值误差按坐标相加。
    lowest=None
    for x,y in product((x0,x1),(y0,y1)):
        bound=endpoint(key,x,y)-cost
        if not bound > TARGET:
            return None
        lowest=bound if lowest is None else lowest.min(bound)
    return lowest.lower()


def verify_region(region):
    key, box = region
    leaves=[]
    todo=[box]
    while todo:
        rect=todo.pop()
        x0,x1,y0,y1=rect
        candidates=sorted(DIRECTIONS,key=lambda k:(2*int(k.split(',')[0])-x0-x1)**2
                                              +(2*int(k.split(',')[1])-y0-y1)**2)
        for candidate in candidates[:16]:
            bound=attempt(candidate,rect)
            if bound is not None:
                leaves.append(bound)
                break
        else:
            x0,x1,y0,y1=rect
            if x1-x0 >= y1-y0:
                mid=(x0+x1)//2
                assert x0<mid<x1, (key,rect)
                todo.extend(((x0,mid,y0,y1),(mid,x1,y0,y1)))
            else:
                mid=(y0+y1)//2
                assert y0<mid<y1, (key,rect)
                todo.extend(((x0,x1,y0,mid),(x0,x1,mid,y1)))
    return min(leaves).str(22),len(leaves)


def verify():
    centers=[tuple(map(int,k.split(','))) for k in DIRECTIONS]
    xs=sorted({x for x,_ in centers})
    ys=sorted({y for _,y in centers})
    assert (xs[0],xs[-1],ys[0],ys[-1])==BOX
    assert set(centers)==set(product(xs,ys))
    def cuts(points):
        return [points[0],*[(a+b)//2 for a,b in zip(points,points[1:])],points[-1]]
    xc,yc=cuts(xs),cuts(ys)
    jobs=[(f'{x},{y}',(xc[i],xc[i+1],yc[j],yc[j+1]))
          for i,x in enumerate(xs) for j,y in enumerate(ys)]
    with ProcessPoolExecutor(max_workers=min(8,cpu_count() or 1)) as pool:
        results=list(pool.map(verify_region,jobs))
    floor=min(arb(v).lower() for v,_ in results)
    assert floor > arb("0.01007")
    print('独立标尺矩形（百万分之一）、有理方向数、叶矩形数、规范差下界',
          BOX,len(DIRECTIONS),sum(n for _,n in results),floor.str(13))


if __name__=='__main__':
    verify()

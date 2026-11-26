"""九枚候选的六时间概率多边形，两个新增点共同份额的精确最大值。"""
from itertools import combinations
from math import prod
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, arb_mat, ctx


def verify(root):
    def read(name):
        return run_path(str(root/('ch39-'+name+'.py')))
    exact=read('six-time-proton-exact-zero')
    previous=read('six-time-positive-support')
    support, pair=exact['SUPPORT'], previous['PAIR']
    heights, p=read('six-time-proton-small-mean')['six_time_support'](root,support)
    ctx.prec=224
    added=[acb.zeta_zero(n).imag for n in pair]
    ctx.prec=320
    scale=read('six-time-height-bands')['scale']
    col=previous['time_column']
    # 六时间核乘公共分母 D 后张成六次以下多项式。用重心插值计算 V^{-1}W，
    # 避免以几乎重合的列直接解未预处理的方程。
    coordinates=[1/(1+(h/scale)**2) for h in (*[heights[n] for n in support],*added)]
    def denominator(y):
        return prod(1+(j*j-1)*y for j in range(2,7))
    def coordinate_coefficient(i,t):
        node=coordinates[i]
        return denominator(node)/denominator(t)*prod((t-u)/(node-u) for j,u in enumerate(coordinates[:7]) if j!=i)
    coeff=arb_mat([[coordinate_coefficient(i,t) for t in coordinates[7:]] for i in range(7)])
    V=arb_mat([col(heights[n],scale) for n in support]).transpose()
    W=arb_mat([col(h,scale) for h in added]).transpose()
    assert all(v.contains(0) for v in (V*coeff-W).entries())
    # 所有余量取非负，末两行分别代表新添的两个概率。
    halfplanes=[(v,-coeff[i,0],-coeff[i,1]) for i,v in enumerate(p)]
    halfplanes += [(arb(0),arb(1),arb(0)),(arb(0),arb(0),arb(1))]
    def intersect(i,j):
        a,b,c=halfplanes[i]
        u,v,w=halfplanes[j]
        det=b*w-c*v
        assert not det.contains(0),(i,j,det,b,c,v,w)
        return (c*u-a*w)/det,(a*v-b*u)/det
    # 除轴线与对角线的公共原点外，所有可行交点即使考虑区间误差也无退化。
    def feasible(x,y,excluded):
        remainder=[a+b*x+c*y for k,(a,b,c) in enumerate(halfplanes) if k not in excluded]
        if any(z<0 for z in remainder):return False
        assert all(z>0 for z in remainder)
        return True
    vertices=[]
    for i,j in combinations(range(9),2):
        x,y=intersect(i,j)
        if feasible(x,y,(i,j)):
            print('多边形顶点，取零行号、两新增概率',i+1,j+1,x.str(17),y.str(17))
            vertices.append((i,j,x,y))
    diagonals=[]
    for i,(a,b,c) in enumerate(halfplanes[:7]):
        assert not (b+c).contains(0)
        e=-a/(b+c)
        if feasible(e,e,(i,)):
            diagonals.append((i,e))
            print('可行正对角交点',i+1,e.str(18))
    assert {(i,j) for i,j,_,_ in vertices}=={(2,6),(2,8),(6,7),(7,8)}
    assert len(diagonals)==1
    row,eta=diagonals[0]
    assert all(eta>x.min(y) for _,_,x,y in vertices)
    a,b,c=halfplanes[row]
    assert row==2 and b<0 and c<0
    assert arb('.22617376849209479562')<eta<arb('.22617376849209479564')
    old=[a+b*eta+c*eta for a,b,c in halfplanes[:7]]
    assert old[row].contains(0) and all(v>0 for j,v in enumerate(old) if j!=row)
    print('最大共同概率和旧支撑概率',eta.str(20),*[v.str(15) for v in old])
    # 利用核对过的统一有理两点分离下界，给所有可行作用的均方下界。
    m=[arb(z)/previous['SCALE'] for z in previous['MIX']]
    floor=eta*arb('.00456315')**2/sum(z*z for z in m)
    assert floor>arb('.00000755845769913')
    print('两点贡献的普适严格下界',floor.str(18))
    # 上一构造有七个精确零，在这张九点字典上的均方完全由新增两点给出。
    shell,lib=read('three-ray-chiral-boundary'),read('mincut-global-winding')
    fixed,angles=read('height-independent-flavor-vacuum'),read('mincut-reference-ckm-action')
    flavor=read('chiral-neutral-flavor-phase')
    ctx.prec=224
    ref_x=arb('-.01046')
    ref_y=ref_x-angles['reference_angles'](flavor)[0]
    tri=lib['terms'](shell)
    def gap(h):
        phases=shell['geo']['phases'](h)
        return (phases[724]/phases[804]).arg()
    gaps=[gap(heights[n]) for n in support]
    locked={q:arb(v) for q,v in exact['FIXED'].items()}
    solved=exact['SOLVED']
    bare={q:lib['coeff'](tri,q,arb(0))[2]/12 for q in (*locked,*solved)}
    def component(q,d):
        return (bare[q]*acb(0,q*d).exp()).real
    columns={q:[arb(1),*fixed['gradient_columns'](lib,tri,q,ref_x,ref_y),
                *(component(q,d) for d in gaps)] for q in bare}
    B=arb_mat([columns[q] for q in solved]).transpose()
    rhs=arb_mat([[int(i==0)-sum(w*columns[q][i] for q,w in locked.items())] for i in range(12)])
    locked.update(zip(solved,B.solve(rhs).entries(),strict=True))
    values=[sum(w*component(q,gap(h)) for q,w in locked.items()) for h in added]
    parts=[eta*v*v for v in values]
    assert arb('.00005107183861498')<sum(parts)<arb('.00005107183861500')
    assert sum(parts)>floor
    print('七零作用的两读数、分项贡献及总均方',*[v.str(17) for v in (*values,*parts,sum(parts))])


if __name__=='__main__':
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)

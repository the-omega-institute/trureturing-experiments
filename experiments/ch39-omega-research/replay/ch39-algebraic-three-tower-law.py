"""由四个钟坐标和径向亏额生成三塔的有限实多项式系统，并核验分支和次数。"""
import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations, product
from math import factorial, prod
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx


# 每一目标依次记录原块、块内升序索引、原谱相对于 c 的符号和目标块。
TAGS = {
    2: ((2,2,1,0),(1,2,1,2),(0,0,-1,1),(0,3,1,0),(2,3,1,2),(0,4,1,1),
        (1,3,1,0),(0,2,-1,0),(0,1,-1,2),(2,4,1,1),(1,1,-1,0),(2,1,-1,2),
        (1,4,1,1),(1,0,-1,2),(2,0,-1,1),(0,5,1,0),(2,5,1,1),(0,6,1,2),(1,5,1,0)),
    4: ((1,2,1,0),(2,2,1,2),(0,0,-1,1),(0,3,1,0),(2,3,1,2),(1,3,1,0),
        (0,4,1,1),(0,2,-1,0),(0,1,-1,2),(1,1,-1,1),(2,4,1,0),(2,1,-1,2),
        (1,4,1,1),(1,0,-1,2),(2,0,-1,1),(0,5,1,0),(2,5,1,1),(0,6,1,0),(1,5,1,2)),
    6: ((0,0,-1,0),(2,2,-1,2),(1,2,1,1),(0,3,1,0),(2,3,1,0),(1,3,1,2),
        (0,4,1,1),(0,2,-1,0),(0,1,-1,2),(2,4,1,1),(1,1,-1,0),(2,1,-1,2),
        (1,4,1,1),(1,0,-1,2),(2,0,-1,1),(0,5,1,0),(2,5,1,1),(0,6,1,2),(1,5,1,0))
}
SECTIONS = {2:(0,2,6,16,10,11), 4:(0,1,7,9,17,12), 6:(0,2,6,9,10,11)}


def geometry(prev):
    outer = sorted(set(permutations((1,5,29))))
    inner = sorted(set(permutations((13,13,23))))
    sites = outer+inner+sorted(v for v in product(range(30),repeat=3)
                if sum(x*x for x in v)==867 and v not in outer+inner)
    rotate = [sites.index((y,z,x)) for x,y,z in sites]
    cycles = [(i,rotate[i],rotate[rotate[i]]) if i!=rotate[i] else (i,)
              for i in range(19) if i==min(i,rotate[i],rotate[rotate[i]])]
    pole = sites.index((17,17,17))
    def dot(a,b):
        return sum(x*y for x,y in zip(a,b))
    edges = {}
    for i,j in combinations(range(19),2):
        x,y = sites[i],sites[j]
        d = 2*(867-dot(x,y))
        if d>804:
            continue
        v = sites[pole]
        a = 867*(867+dot(v,x)+dot(x,y)+dot(y,v))
        b = -17*sum(v[t]*(x[(t+1)%3]*y[(t+2)%3]-x[(t+2)%3]*y[(t+1)%3]) for t in range(3))
        edges[i,j] = (d,a,b)
    groups = prev["orbits"]({"sites":sites,"edges":[(i,j,d,None) for (i,j),(d,_,_) in edges.items()]})
    return sites, cycles, edges, groups, pole


def polynomial_system(root, charge, inputs=None):
    """返回五个输入、实未知量、零多项式、严格正多项式和按轨道排序的 45 个输出。"""
    import sympy as s
    prev=run_path(str(root/"ch39-three-even-positive-towers.py"))
    _,cycles,edges,groups,pole=geometry(prev)
    coords = s.symbols("cx sx cy sy delta", real=True) if inputs is None else tuple(map(s.sympify,inputs))
    cx,sx,cy,sy,delta=coords
    native=[1-delta if any(pole in e and set(e)&{0,3,4} for e in g) else 1 for g in groups]
    circle=[cx*cx+sx*sx-1,cy*cy+sy*sy-1]
    residue=charge%14
    q=min(residue,14-residue)
    if q==0 or q%2:
        return {"inputs":coords,"unknowns":(),"zero":circle,"positive":(),"weights":native}
    ph={724:cx+s.I*sx,804:cy+s.I*sy}
    num,variable,_=prev["DATA"][q]
    w=s.symbols("w0:19", real=True)
    base=[s.Rational(v,10**15) for v in num]
    fitted=base.copy()
    for i,a in zip(variable,w):
        fitted[i]=a
    unit={}
    for (i,j),(d,a,b) in edges.items():
        # 偶荷消去单位化的正平方根，d=200 的负号也消失。
        squared=(a+s.I*s.sqrt(3)*b)**2/(a*a+3*b*b)
        unit[i,j]=s.S.One if pole in (i,j) else squared**(q//2)*ph.get(d,1)**q
    omega=(-1+s.I*s.sqrt(3))/2
    def blocks(ws):
        L=s.zeros(19)
        for weight,g in zip(ws,groups):
            for i,j in g:
                L[i,i]+=weight
                L[j,j]+=weight
                L[i,j]-=weight*unit[i,j]
                L[j,i]-=weight*s.conjugate(unit[i,j])
        result=[]
        for k in range(3):
            use=[g for g in cycles if k==0 or len(g)==3]
            result.append(s.ImmutableMatrix([[sum(L[row[0],j]*omega**((k*t)%3)
                         for t,j in enumerate(col)) for col in use] for row in use]))
        return result
    native_blocks,fit_blocks=blocks(native),blocks(fitted)
    roots=[s.symbols(f"u{k}_0:{B.rows}",real=True) for k,B in enumerate(native_blocks)]
    targets=s.symbols("y0:19",real=True)
    def det_at(B,z):
        return s.re(s.Determinant(z*s.eye(B.rows)-B),evaluate=False)
    eq=circle
    pos=[]
    for B,us in zip(native_blocks,roots):
        eq.extend(det_at(B,u) for u in us)
        pos.extend(v-u for u,v in zip(us,us[1:]))
    c=15-delta/3
    for i,(k,j,sign,dest) in enumerate(TAGS[q]):
        u,y=roots[k][j],targets[i]
        eps=-1 if i<2 else 1
        factor=1 if i<3 else 32
        eq.extend((y*(y-c)-eps*factor*sign*u*(u-c),det_at(fit_blocks[dest],y)))
        pos.extend((sign*(u-c), (c-2*y) if i==0 else (2*y-c)))
    pos.extend(v-u for u,v in zip(targets,targets[1:]))
    for i in variable:
        half=s.Rational(2,10**7)
        pos.extend((half+(fitted[i]-base[i]),half-(fitted[i]-base[i])))
    return {"inputs":coords, "unknowns":tuple(u for row in roots for u in row)+targets+w,
            "zero":eq, "positive":pos, "weights":fitted}


def rational_reference():
    denominator = 1_044_988
    inputs = []
    for numerator in (-5_469, 72_843):
        norm = denominator**2 + numerator**2
        inputs.extend((Fraction(denominator**2 - numerator**2, norm),
                       Fraction(2 * denominator * numerator, norm)))
    inputs.append(Fraction(1, 2**98))
    return inputs


def mixed_bound():
    coefficients=[1]
    for degree,cap in [(j,b) for j in range(1,7) for b in (j,j-1,j-1)]+[(7,6)]:
        coefficients=[a+b for a,b in zip([cap*v for v in coefficients]+[0],
                                         [0]+[(degree-cap)*v for v in coefficients])]
    assert sum(coefficients)==factorial(7)*factorial(6)**2
    return sum(coefficients[:11])


def verify(root):
    prev=run_path(str(root/"ch39-three-even-positive-towers.py"))
    old=run_path(str(root/"ch39-three-tower-phase-continuation.py"))
    old["verify"](root,root/"ch39-three-even-positive-towers.py")
    shell=run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec=288
    sites,cycles,edges,groups,pole=geometry(prev)
    assert sites==shell["geo"]["sites"] and pole==shell["soft"]["pole"]
    assert set(shell["soft"]["star"]["orbit"])=={0,3,4}
    assert groups==prev["orbits"](shell["geo"])
    bases=old["sectors"](sites)
    h=acb.zeta_zero(prev["NUMBER"]).imag
    clocks=shell["geo"]["phases"](h)
    delta=arb(23)/25*arb(2)**-20/(h*h+shell["soft"]["point"]["study"]["offset"])
    c=15-delta/3
    bound=mixed_bound()
    assert bound==2612295360
    for q in TAGS:
        num,var,perm=prev["DATA"][q]
        original=shell["laplacian"](q,h)
        phases={}
        for (i,j),(d,a,b) in edges.items():
            u=acb(a,arb(3).sqrt()*b)
            z=(u*u/(a*a+3*b*b))**(q//2)*clocks.get(d,1)**q
            phases[i,j]=acb(1) if pole in (i,j) else z
            assert (phases[i,j]+original[i,j]/abs(original[i,j])).contains(0)
        weights=[arb(n)/10**15 for n in num]
        geom={e:v[0] for e,v in edges.items()}
        fitted,axes,forms=old["data"](groups,geom,phases,weights,q)
        def rows(matrix):
            out=[]
            for k,T in enumerate(bases):
                block=T.conjugate().transpose()*matrix*T
                values=sorted((z.real for z in block.eig(multiple=True)),key=float)
                out.extend((v,k,i) for i,v in enumerate(values))
            return sorted(out,key=lambda row:float(row[0]))
        native,fit=rows(original),rows(fitted)
        dn=18*q*arb("1e-9")+4*arb("1e-12")
        df=sum(old["rownorm"](v) for v in axes)*arb("1e-9")
        inc=[0]*19
        for g in var:
            for i,j in groups[g]:
                inc[i]+=1
                inc[j]+=1
        change=df+2*max(inc)*arb("2e-7")
        assert all(v[0]-u[0]>2*change for u,v in zip(fit,fit[1:]))
        marks=sorted(native,key=lambda u:float(u[0]*abs(u[0]-c)))
        cw=c+arb(0,arb("1e-12")/3)
        intervals=[]
        for v,_,_ in marks:
            z=v+arb(0,dn)
            assert z<cw or z>cw
            intervals.append(z*abs(z-cw))
        assert all(a<b for a,b in zip(intervals,intervals[1:]))
        observed=[]
        for out,i in enumerate((*perm,*range(3,19))):
            z,k,j=marks[i]
            observed.append((k,j,1 if z>c else -1,fit[out][1]))
        assert tuple(observed)==TAGS[q]
        ys=old["targets"](intervals,cw,perm)[0]
        assert all(a<b for a,b in zip(ys,ys[1:]))
        # 在有理隔离盒直接核验导数单射。原开盒证明给出的解位于此盒内。
        eigrows=old["eigensystem"](fitted,bases)
        jac=arb_mat([[old["rayleigh"](u,forms[g]) for g in var] for _,u,_ in eigrows])
        inv=jac.inv()
        err=arb_mat([[12*change/(gap-2*change)+
                      (3*q*arb("1e-9") if any(geom[e] in (724,804) for e in groups[g]) else 0)
                      for g in var] for _,_,gap in eigrows])
        invabs=arb_mat([[abs(inv[i,j]) for j in range(19)] for i in range(19)])
        contraction=old["rownorm"](invabs*err)
        assert contraction<1
        section=set(SECTIONS[q])
        assert all(len(section.intersection(cyc))==int(len(cyc)==3) for cyc in cycles)
        selected=[g for g in var if any(set(e)<=section for e in groups[g])]
        assert len(selected)==9 and all(sum(set(e)<=section for e in groups[g])==1 for g in selected)
        tags=Counter((k,sign,dest, i if i<3 else 3) for i,(k,_,sign,dest) in enumerate(TAGS[q]))
        partitions=factorial(7)*factorial(6)**2//prod(factorial(v) for v in tags.values())
        radicals=16+sum(t[2]!=s for t,s in zip(TAGS[q][:3],(-1,-1,1)))
        field=partitions*2**radicals
        print("荷、隔离盒压缩率、九轨道、原根分组数、平方根数、目标域次数界、总次数界",
              q,contraction.str(12),selected,partitions,radicals,field,bound*field)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root",type=Path,nargs="?",default=Path(__file__).resolve().parent)
    parser.add_argument("--system",type=int,help="输出指定整数荷的符号实多项式系统")
    parser.add_argument("--rational-reference", action="store_true",
                        help="直接使用三塔角盒内的整数半角及有理亏额生成方程")
    args=parser.parse_args()
    if args.rational_reference and args.system is None:
        parser.error("--rational-reference 需与 --system 指定的整数荷合用")
    if args.system is None:
        verify(args.root)
    else:
        inputs = rational_reference() if args.rational_reference else None
        system=polynomial_system(args.root, args.system, inputs)
        for label,terms in system.items():
            print(label)
            for term in terms:
                print(term)

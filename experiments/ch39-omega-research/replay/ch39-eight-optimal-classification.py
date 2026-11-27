"""列举能同时消去七壳弱差的全部最小八站支撑。"""
from itertools import combinations, product
from pathlib import Path
from runpy import run_path
from flint import acb, acb_mat, arb, ctx


def adj(a):
    return a.conjugate().transpose()


def matrix_rows(a, rows):
    return acb_mat([a.tolist()[i] for i in rows])


def exclusion_test(base, gram, scale, inv, drift, rad):
    """利用六个壳行之间的相位关系给出局部 Gram 必要条件。"""
    w = acb_mat([[s*z for z in row] for s,row in zip(scale,inv.tolist(),strict=True)])
    y = base+acb_mat([[s*z for z in row] for s,row in zip(scale,drift.tolist(),strict=True)])
    d = gram-adj(y)*y
    a,b,c=d[0,0].real,d[0,1],d[1,1].real
    lo=(a+c-((a-c)**2+4*abs(b)**2).sqrt())/2
    # 任意测试向量均可。用区间中点的有理近似选方向，依其实际范数归一化。
    v=acb_mat([[-b.mid()],[(a-lo).mid()]])
    length=sum(abs(z)**2 for z in v.entries()).sqrt()
    if not length>0:
        v=acb_mat([[1],[0]])
    else:
        v=v/length
    coeff=(adj(y*v)*w).entries()
    err=2*sum(abs(z)*r for z,r in zip(coeff,rad,strict=True))
    result=(-(adj(v)*d*v)[0,0].real-err).lower()
    assert result.is_finite()
    return result


def experiment(root, fifth, columns, levels, *, observe=None):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    ctx.prec=256
    graph=lib("three-ray-chiral-boundary")
    native=lib("native-rh-ray-pvm")
    gram_native=native["gram_cholesky"](root)
    frame=lib("cyclic-two-characters-nine")["characters"](graph["geo"]["sites"])[0][0]
    orbits=[[i for i in range(19) if abs(frame[i,j])>arb(".5")] for j in range(7)]
    assert orbits[6]==[14]
    f,b,t=levels
    x=fifth*columns*acb_mat([[(b-f).sqrt(),0],[0,(t-f).sqrt()]])
    source=adj(frame)*x
    rows,weights=[],[]
    for number in range(2_306_888_100_092,2_306_888_100_099):
        _,basis=lib("adjacent-pvm-transport")["spectrum"](graph,6,number)
        rows.append(adj(basis))
        weights.append([z.real for z in native["readout"](basis,range(19),[14,8,7,6],gram_native[3],3)])
    first=[(row*frame).tolist()[0] for row in rows]
    full=acb_mat(first).inv()
    planck,mz=arb("1.22089e19"),arb("91.1876")
    limits=(700,700,2000,2000,700,700,700)
    rad=[(arb(q)/planck-f).sqrt() for q in limits]
    # 六轨支撑迫使每壳只有第一条候选能低于 3000 GeV，细截断界由弱核给出。
    charts=[]
    for drop in [2,3,1,5,0,4,6]:
        ix=[j for j in range(7) if j!=drop]
        inv=acb_mat([first[j][:6] for j in ix]).inv()
        drift=(-inv*acb_mat([[first[j][6]] for j in ix]))*matrix_rows(source,[6])
        charts.append((inv,drift,[rad[j] for j in ix]))
    centers=list(map(arb,(*lib("six-cyclic-weak-horizon")["Q"],"506.579672393")))
    width=arb(".00005")
    radii=[((q+arb(0,width.upper()))/planck-f).sqrt() for q in centers]
    pole=matrix_rows(source,[6])
    length=sum(abs(z)**2 for z in pole.entries()).sqrt()
    direction=pole/length
    sides=[abs(full[6,j])*r for j,r in enumerate(radii)]
    # 第八、九壳组成一组，其余五壳为另一组。全部 32 个成功支撑共享这一相位分组。
    mask=6
    left=sum(sides[1:3])
    right=sum(z for j,z in enumerate(sides) if j not in (1,2))
    assert abs(left-right)<length<left+right
    xp=(length**2+left**2-right**2)/(2*length)
    yp=(left**2-xp**2).sqrt()
    phases=[acb(xp,yp)/left if j in (1,2) else acb(length-xp,-yp)/right for j in range(7)]
    projected_target=acb_mat([[r*p*abs(full[6,j])/full[6,j]]
                              for j,(r,p) in enumerate(zip(radii,phases,strict=True))])
    target=full*projected_target*direction
    reference=arb("-.074256577194663365546").union(arb("-.074256577194663365544"))
    p0,p1=(mz/planck)**2,(arb("1.176397e14")/planck)**2
    weak=native["weak_series"](56)
    large=lib("full-neutral-free-cutoff")["weak_column"]
    faces=[(reference+w[0]*(weak(((q-width)/mz)**2)-reference),
            reference+w[0]*(weak(((q+width)/mz)**2)-reference)) for q,w in zip(centers,weights,strict=True)]
    count=0
    classified=[]
    rotated_gaps={}
    least_gap,least_margin=arb(10),arb(10)
    for sizes in product((1,2,3),repeat=6):
        if sum(sizes)!=8:
            continue
        scale=[(arb(3)/s).sqrt() for s in sizes]
        for groups in product(*(combinations(orbit,s) for orbit,s in zip(orbits[:6],sizes,strict=True))):
            count+=1
            sites=[i for g in groups for i in g]
            local=matrix_rows(x,sites)
            gram=adj(local)*local
            base=acb_mat([[sum(x[i,h] for i in g)/arb(len(g)).sqrt()-a*source[k,h] for h in (0,1)]
                          for k,(g,a) in enumerate(zip(groups,scale,strict=True))])
            for inv,drift,rad6 in charts:
                gap=exclusion_test(base,gram,scale,inv,drift,rad6)
                if gap>0:
                    least_gap=min(least_gap,gap)
                    key=tuple(sorted(sites))
                    if key in {(0,2,4,5,6,12,15,17),(0,1,2,3,7,13,16,18)}:
                        rotated_gaps[key]=gap
                    break
            else:
                avg=base+acb_mat([[scale[k]*target[k,h] for h in (0,1)] for k in range(6)])
                residual=gram-adj(avg)*avg
                assert residual[0,0].real>0,groups
                a=residual[0,0].real.sqrt()
                b=residual[0,1]/a
                d2=residual[1,1].real-abs(b)**2
                assert d2>0,groups
                completion=acb_mat([[a,b],[0,d2.sqrt()]])
                image=acb_mat(x)
                extra=0
                for k,g in enumerate(groups):
                    for i in g:
                        for h in (0,1):
                            image[i,h]=avg[k,h]/arb(len(g)).sqrt()
                    # t 个 1 和一个 -t 的实向量形成该轨道的零和标准正交组。
                    for t in range(1,len(g)):
                        scale0=arb(t*(t+1)).sqrt()
                        for j,i in enumerate(g[:t+1]):
                            for h in (0,1):
                                image[i,h]+=(1 if j<t else -t)*completion[extra,h]/scale0
                        extra+=1
                assert extra==2
                margin=arb(1)
                for row,w,face in zip(rows,weights,faces,strict=True):
                    hi=[(f+sum(max(arb(0),abs(z).lower())**2 for z in r)).lower()
                        for r in (row*image).tolist()[1:]]
                    assert all(z.is_finite() for z in hi)
                    high=min(hi)
                    assert high**2>p0,groups
                    err=(large(high**2,p0,p1)-reference).abs_upper()*sum(abs(z).upper() for z in w[1:])
                    m=min((face[0]-err).lower(),(-face[1]-err).lower())
                    assert m.is_finite()
                    margin=min(margin,m)
                assert margin>0,groups
                least_margin=min(least_margin,margin)
                classified.append((sizes,tuple(sorted(sites)),mask))
                if observe is not None:
                    observe(tuple(sorted(sites)), image)
    from collections import Counter
    counts=Counter(s for s,sites,m in classified)
    assert {mask for _,_,mask in classified}=={6}
    assert count==12393 and counts=={(2,2,1,1,1,1):27,(3,1,1,1,1,1):5},counts
    # 固定同一个第五壳源、第三枢轴和弱读数，只在坐标站集上施加 (a,b,c) ↦ (b,c,a)。
    geom_sites=graph["geo"]["sites"]
    rot=[geom_sites.index((b,c,a)) for a,b,c in geom_sites]
    def rotate(s):
        return tuple(sorted(rot[i] for i in s))
    winners={s for _,s,_ in classified}
    representatives={min(s,rotate(s),rotate(rotate(s))) for s in winners}
    occupancy=Counter()
    for s in representatives:
        orbit={s,rotate(s),rotate(rotate(s))}
        assert len(orbit)==3 and rotate(rotate(rotate(s)))==s
        occupancy[len(winners & orbit)]+=1
    assert occupancy=={1:8,2:12}
    example=(1,3,4,5,8,9,10,11)
    turned=rotate(example)
    third=rotate(turned)
    assert turned==(0,2,4,5,6,12,15,17) and third==(0,1,2,3,7,13,16,18)
    assert example in winners and turned not in winners and third not in winners
    print("支撑循环轨道内恰有一、二、三个可行支撑的轨道数",occupancy[1],occupancy[2],occupancy[3],flush=True)
    assert set(rotated_gaps)=={turned,third}
    print("可行八站支撑及两次被排除的坐标循环",example,turned,third,flush=True)
    print("两次坐标循环的方向 Gram 排除余量",rotated_gaps[turned].str(12),rotated_gaps[third].str(12),flush=True)
    supports = [set(sites) for _,sites,_ in classified]
    assert len({tuple(sorted(s)) for s in supports})==32
    assert set.intersection(*supports)==set()
    assert set.union(*supports)==set(range(19))-{14}
    print("保极候选、可实现分布、排除最小余量、构造最小面余量",count,counts,least_gap.str(12),least_margin.str(12),flush=True)
    for _,sites,mask in classified:
        print("八站支撑、相位组掩码",sites,mask)
    return winners

    return winners


if __name__=="__main__":
    root=Path(__file__).resolve().parent
    source=run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    experiment(root,*source["verify"](root,return_operator="spectral"))

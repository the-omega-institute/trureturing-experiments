"""验证保持坐标三循环的正边电导逐层产生实际零点的完整荷二十六重层目标。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, arb_mat, ctx


NUMBER=2_306_888_100_085
CENTERS=(
    6511705252361918, 2825066589748832, 3097598382619734, 1811590548198980,
    2236427605769425, 1568977481681329, 4245867285407642, 1920281309411106,
    6288885771809051, 3308468971584334, 2696792337624010, 3022811559624942,
    6147868866236653, 4752394317835150, 5657679796329590, 2810683110530149,
    3290302238718345, 3069152627257324, 4383820938881652, 2507267652964743,
    1997189925178908, 2764138575404864, 2710079356092523, 5245609731754697,
    4762998577074440, 2132061974439302, 700509728910781, 2529051145816817,
    767754320155027, 1772391359788810, 2783510464520065, 4763813294699755,
    3659603897414109, 1802224449312906, 2620097319394101, 1998037045889745,
    2808351834950756, 1692477810597548, 5005228292225302, 2367672678708062,
    1754120353612798, 2384754715600978, 1379134482410167, 2682023107647060,
    185728620783454,
)
VARIABLE=(8, 35, 19, 6, 11, 28, 37, 0, 43, 10, 40, 21, 15, 39, 20, 30)

ALL_CENTERS=(
    8823140598458227, 1272446501112676, 1170800482052140, 1123461284052940,
    6196797338200720, 311618678040991, 9889602805145154, 5341784489335363,
    240411529682013, 528023277808278, 2749733033154498, 5403081118092979,
    3342978306546984, 692831631979121, 5128909682871235, 7829288801736639,
    1875261112912928, 2775913824651806, 609083476265640, 3816204310836278,
    4061567904163028, 3381961659714570, 1791896225166656, 1422250191051521,
    4733588246832641, 4123293156692201, 3951049773674640, 2013493015362617,
    2954393943688120, 1328559907980052, 2475659170810364, 4530657739265957,
    791235795319092, 2531934248774394, 6293570866927919, 1347470220437038,
    4729942498066868, 1191539494153184, 601013138709505, 5239148556415989,
    854166504956352, 129954691618509, 631749323000149, 173837895049385,
    1384402190232944,
)
ALL_VARIABLE=(19, 6, 41, 35, 34, 14, 8, 43, 40, 11, 42, 13, 10, 15, 17, 5, 9, 33, 12)


def verify(root):
    shell=run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec=288
    h=acb.zeta_zero(NUMBER).imag
    native=shell["laplacian"](2,h)
    neutral=shell["laplacian"](0,h)
    c=sum(neutral[j,j].real for j in range(19))/18
    sites=shell["geo"]["sites"]
    edges={(i,j) for i,j,_,_ in shell["geo"]["edges"]}
    p=[sites.index((y,z,x)) for x,y,z in sites]
    groups=[]
    phases={}
    contrasts={}
    unseen=set(edges)
    while unseen:
        edge=min(unseen)
        orbit=[]
        for _ in range(3):
            i,j=edge
            orbit.append(edge)
            amplitude=abs(native[i,j])
            phases[edge]=-native[i,j]/amplitude
            contrasts[edge]=amplitude
            edge=tuple(sorted((p[i],p[j])))
        assert len(set(orbit))==3 and edge==orbit[0]
        unseen.difference_update(orbit)
        groups.append(orbit)
    original=native.eig(multiple=True)
    assert len(original)==19 and all(z.imag.contains(0) and z.real>0 for z in original)
    folded=sorted((z.real*abs(z.real-c) for z in original),key=float)
    assert all(a<b for a,b in zip(folded,folded[1:])) and 32*folded[3]>c*c/4
    target=[(c+(c*c+128*z).sqrt())/2 for z in folded[3:]]
    all_target=[(c-(c*c-4*folded[1]).sqrt())/2,
                (c+(c*c-4*folded[2]).sqrt())/2,
                (c+(c*c+4*folded[0]).sqrt())/2,*target]
    assert all(v<w for v,w in zip(all_target,all_target[1:]))
    assert all(m<9 for m in folded[:3])
    for name,center_values,variables,goals in (
        ("只固定十六个重层",CENTERS,VARIABLE,target),
        ("同时精确保留原三个轻层",ALL_CENTERS,ALL_VARIABLE,all_target)):
        print(name)
        certify(center_values,variables,goals,groups,phases,contrasts,c,folded)


def certify(center_values,variables,goals,groups,phases,contrasts,c,folded):
    n=len(goals)
    offset=19-n
    assert len(groups)==len(center_values)==45 and len(set(variables))==n
    centers=[arb(v)/10**15 for v in center_values]
    base=acb_mat(19,19)
    for w,orbit in zip(centers,groups):
        for i,j in orbit:
            entry=-w*phases[i,j]
            base[i,i]+=w
            base[j,j]+=w
            base[i,j],base[j,i]=entry,entry.conjugate()
    values,vec=base.eig(right=True,algorithm="rump")
    order=sorted(range(19),key=lambda k:float(values[k].real))
    eig=[values[k].real for k in order]
    assert all(z.imag.contains(0) and z.real>0 for z in values)
    gap=min((b-a).lower() for a,b in zip(eig,eig[1:]))
    assert gap>arb(".01")
    derivative=arb_mat(n,n)
    for mode,k in enumerate(order[offset:]):
        norm=sum(abs(vec[i,k])**2 for i in range(19))
        assert norm>0
        for j,g in enumerate(variables):
            derivative[mode,j]=sum(abs(vec[i,k]-phases[i,l]*vec[l,k])**2 for i,l in groups[g])/norm
    pre=derivative.inv()
    invnorm=max(sum(abs(pre[i,j]) for j in range(n)).upper() for i in range(n))
    residual=arb_mat([[a-b] for a,b in zip(eig[offset:],goals)])
    correction=max(abs(v).upper() for v in (pre*residual).entries())
    radius=arb("1e-10")
    # 每个三边轨道的算子范数 ≤6；全盒所有特征值仍为单根。
    matrix_change=6*n*radius
    min_gap=gap-2*matrix_change
    # 沿任一直线投影导数的迹范数 ≤2 ||L'||/gap。
    entry_error=12*matrix_change/min_gap
    contraction=invnorm*n*entry_error
    assert contraction<1 and correction+contraction*radius<radius
    print("中点最小谱隙、逆雅可比无穷范数、预乘中点残差、整盒压缩率上界",
          *[v.str(15) for v in (gap,invnorm,correction,contraction)])
    first=[v+arb(0,matrix_change) for v in eig[:3]]
    lowmass=[v*abs(v-c) for v in first]
    assert all(m<32*folded[3] for m in lowmass)
    if n==19:
        assert all(z<9 for z in lowmass)
    else:
        assert lowmass[0]<9 and all(z>10 for z in lowmass[1:])
    print("前三个本征值、在精确解处前三个未排序的折叠质量",
          *[z.str(12) for z in (*first,*lowmass)])
    limit=arb(0)
    for g,(w,orbit) in enumerate(zip(centers,groups)):
        unknown=w+arb(0,radius) if g in variables else w
        assert unknown>0
        for e in orbit:
            ref=contrasts[e]
            limit=limit.max((unknown/ref).max(ref/unknown))
    assert limit<(arb("6.511706") if n==16 else arb("9.889604"))
    print("逐边双向反差与对数的可达上界",limit.upper().str(16),limit.log().upper().str(16))
    print("精确匹配的第4、第19个折叠质量目标",(32*folded[3]).str(16),(32*folded[-1]).str(16))


if __name__=="__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)

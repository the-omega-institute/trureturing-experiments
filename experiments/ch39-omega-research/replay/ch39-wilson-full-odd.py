"""指定三座偶荷逆谱塔时，物质奇荷一阶迹的最小代价。"""
from json import loads
from pathlib import Path
from runpy import run_path
import sys
from flint import arb,acb,arb_mat,ctx


def sub(A,rows,cols):
    return arb_mat([[A[i,j] for j in cols] for i in rows])


def geometry(root,helper,record):
    old=run_path(str(helper))
    B,W,t0,_,_,edges,groups,terms,H,ed,base,shell=old['geometry'](root)
    # 利用整数入射关系核查 177 个坐标覆盖全部允许的平均响应。
    old['certify_coordinates'](root,record,B,edges,groups,terms,ed,shell)
    tomography=run_path(str(root/'ch39-mincut-real-triangle-tomography.py'))
    labels=tomography['enumerate_faces']()
    towers=run_path(str(root/'ch39-three-even-positive-towers.py'))
    towers['verify'](root)
    ctx.prec=512
    # 直接用整数面标记及高精度双钟重建矩阵，避免消元时放大中间舍入。
    R=tomography['R']
    def theta(d):
        return (1-arb(d)/(2*R)).acos()
    h=acb.zeta_zero(2306888100085).imag
    angles=[(theta(d)-theta(244))*h for d in (724,804)]
    precise=arb_mat(270,179)
    for j,g in enumerate(groups):
        idx=terms[g[0]][0]
        m,(a,b)=labels[tuple(sorted({v for i in idx for v in ed[i]}))]
        native=acb(a,arb(R).sqrt()*b)
        loop=native/abs(native)*acb(0,m[0]*angles[0]+m[1]*angles[1]).exp()
        powers=[1-(loop**q).real for q in range(1,7)]
        for i,e in enumerate(edges):
            count=sum(e[0] in terms[f][0] for f in g)
            for q,v in enumerate(powers):
                precise[45*q+i,j]=count*v
    assert all(x.overlaps(y) for x,y in zip(B.entries(),precise.entries()))
    C,rows=record['C'],record['R']
    U=sub(precise,range(270),C)*sub(precise,rows,C).inv()
    # 只对逆谱证书中的活动变量加盒宽，其他电导取固定有理值。
    orbit=towers['orbits'](shell['geo'])
    mapped={e:j for j,g in enumerate(orbit) for e in g}
    goal=[]
    for q in (2,4,6):
        weights,active,_=towers['DATA'][q]
        for g in edges:
            j=mapped[ed[g[0]]]
            assert all(mapped[ed[i]]==j for i in g)
            v=arb(weights[j])/10**15+(arb(0,arb('2e-11')) if j in active else 0)
            goal.append([v/base[ed[g[0]]]-1])
    target=arb_mat(goal)
    even=[45*q+i for q in (1,3,5) for i in range(45)]
    odd=[45*q+i for q in (0,2,4) for i in range(45)]
    return U,old['trace_rows'](U,W,t0),sub(U,even,range(177)),sub(U,odd,range(177)),target,W,t0


def verify(root,helper):
    data=loads(Path(__file__).with_suffix('.json').read_text())
    U,L,E,O,target,W,t0=geometry(root,helper,data)
    even=[1+sum(W[i]*target[45*q+i,0] for i in range(45))/t0 for q in range(3)]
    assert max(t.upper() for t in even)<4
    assert min((1+t).lower() for t in target.entries())>arb('.125')
    # 假定两条物质奇荷迹比均不超过八百万，第一组对偶乘子保证第五荷也有界。
    K=sub(L,[0,2,1,3,5],range(177))
    y=arb_mat([[arb(z)/10**12 for z in data['y']]])
    mu=arb_mat(1,270)
    for i,z in data['mu']:
        assert z>0
        mu[0,i]=arb(z)/10**12
    error=-sub(L,[4],range(177))-y*K-mu*U
    eps=sum(abs(x).upper() for x in error.entries())
    h=arb(8000000)
    trace_bounds=[arb(0).union(h),arb(0).union(h),*even]
    upper=1-sum(y[0,i]*(z-1) for i,z in enumerate(trace_bounds))+sum(mu.entries())
    # 正性给出 ||x||∞ ≤ 46 max(h,T5/T0)；再代入此式以消去 T5。
    assert t0/min(z.lower() for z in W)<46 and eps*46<arb('.000001')
    assert upper.upper()+46*eps*h<10**8 and upper.upper()/(1-46*eps)<10**8
    bound=46*10**8
    kin=data['kin']
    free=sorted(set(range(177))-set(kin))
    assert len(kin)==len(set(kin))==135
    P=sub(O,range(135),kin)*sub(E,range(135),kin).inv()
    G=sub(O,range(135),free)-P*sub(E,range(135),free)
    f=P*target+arb_mat([[1] for _ in range(135)])
    weight=arb_mat([[w/t0 for w in W]])
    k=arb_mat([(weight*sub(G,range(45*q,45*q+45),range(42))).entries() for q in range(3)])
    b=[(weight*sub(f,range(45*q,45*q+45),[0]))[0,0] for q in range(3)]
    # 真正的三个奇荷边比为 f+Gz，z 是 42 个实际主元增量。
    m=arb_mat(1,135)
    for i,z in data['lower']['mu']:
        assert z>=0
        m[0,i]=arb(z)/10**15
    a=arb(data['lower']['lam'])/10**15
    assert 0<=a<=1
    slope=sub(k,[0],range(42))*a+sub(k,[1],range(42))*(1-a)
    norm=sum(abs(v).upper() for v in (slope-m*G).entries())
    lb=a*b[0]+(1-a)*b[1]-(m*f)[0,0]-bound*norm
    near=arb_mat([[arb(v)/10**7] for v in data["close"]])
    near_floor=min(v.lower() for v in (f+G*near).entries())
    near_trace=[b[i]+(sub(k,[i],range(42))*near)[0,0] for i in range(3)]
    print("逼近下确界的见证最小奇荷边比与三种奇荷迹比",near_floor.str(14),*[v.str(17) for v in near_trace],flush=True)
    assert near_floor>arb(".00114")
    x=arb_mat([[arb(v)/10**5] for v in data['primal']])
    odd=f+G*x
    floor=min(v.lower() for v in odd.entries())
    actual=[b[i]+(sub(k,[i],range(42))*x)[0,0] for i in range(3)]
    print("归一化极大迹的严格下界",lb.lower().str(20),"两组对偶残差",eps.str(8),norm.str(8),flush=True)
    print("正见证的最小奇荷边比与荷一、三、五迹比",floor.str(16),*[t.str(20) for t in actual],flush=True)
    assert arb("7729313.16")<lb.lower()<max(t.upper() for t in near_trace[:2])<arb("7729313.84")
    # 在这两种物质奇荷里至少一个的最大原始本征值超过迹的十九分之一。
    level = lb.lower()*t0/19
    assert t0 > arb("269.999") and t0/18 < 15 and level > arb("1.09837e8")
    mass = level*(level-t0/18)*arb("2.435e18")/724**2
    assert mass > arb("5.60e28")
    assert max(t.upper() for t in actual[:2])<arb("7729329.61")
    assert floor>arb(".03614")


if __name__=='__main__':
    here=Path(__file__).resolve().parent
    root=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else here
    del sys.argv[1:]
    verify(root,root/'ch39-wilson-odd-preservation.py')

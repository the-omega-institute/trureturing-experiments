"""固定联合见证的所有正整数中心偶倍率与质量分母，覆盖全部允许的交点能量。"""
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx


def dot(a,b):
    return sum(x*y for x,y in zip(a,b,strict=True))


def verify(root):
    physical = run_path(str(root/'ch39-first-three-mode-observables.py'))
    model, initial = physical['model_inputs']()
    ctx.prec = 128
    joint = physical['joint']
    chiral = joint['chiral']
    height = acb.zeta_zero(2_306_888_100_085).imag
    rows = chiral['source']['spectral_supplier'](model)(height)
    layers = {q:chiral['ordered']([h for h,z in rows if z==q]) for q in (0,*model['charges'])}
    verify_scale_layers(model, initial, joint, layers, arb('31.98396690'), arb('32.01596692'), neighbors=True)


def verify_scale_layers(model, initial, joint, layers, kmin, kmax, *, neighbors=False):
    gauge = joint['gauge']
    sm = gauge['sm']
    unit=arb('2.435e18')/724**2
    # 上一证书的全域热核界给等效倍率 k (724/d)^2 的必要外区间。
    a=(arb('246.22')**2/(arb('1.22089e19')*32*unit*layers[2][15]))**2
    c=(9*arb(3).log()-12*arb(2).log())/2
    nominal_curve=c*a*a/(arb.pi()*arb('1e-122'))
    # 在任一合格牛顿根 x 处，E2<=1 给出不依赖两个尺度的 x 下界。
    # 精确配平的场内容固定权重，移动各层质量只改变热核参数。
    vacuum=joint['vacuum']
    template=vacuum['minimal']['fields'](layers)
    weights=[w+(arb(-g)/6 if g<0 else arb(g)/12) for g,w,z in template]
    weights.extend(w for w,z in vacuum['low_sector'](4,a))
    xmin=arb('.006')
    assert 2*arb.pi()/sum(w for w in weights if w>0)>xmin
    assert nominal_curve/4 < arb('.998')
    assert 4*nominal_curve*(1-2*a/(xmin*c)) > arb('1.002')
    assert 2*a/(xmin*c)<arb('1e-57')
    assert kmin < 32*(nominal_curve*(1-arb('1e-57'))/arb('1.002')).root(4)
    assert 32*(nominal_curve/arb('.998')).root(4)<kmax
    assert kmax*(arb(127)/724)**2<1
    matching=[]
    for d in range(702,757):
        candidate=(kmin.union(kmax)*(arb(d)/724)**2).unique_fmpz()
        if candidate is not None:
            matching.append((int(candidate),d))
        assert (kmax-kmin)*(arb(d)/724)**2 < 1
    assert matching==[(32,724)],matching
    even=[*layers[0],*(h for q in (2,4,6) for h in layers[q][3:])]
    limit=(min(v.upper() for v in even)*unit*kmax/model['mz']).log()
    assert limit<arb('29.4')
    groups=(1,3)
    times={q:[(unit*v/model['mz']).log() for v in layers[q][3:]] for q in groups}
    first=min(v.lower() for a in times.values() for v in a)
    last=max(v.upper() for a in times.values() for v in a)
    assert first<arb('27.013') and last<30
    assert all(a<b for ts in times.values() for a,b in zip(ts,ts[1:]))
    example=gauge['two']['thresholds'](model,{1:[0,0,0,arb(1)],3:[0,0,0,arb(2)]},False)
    one=[r[1] for r in example];two=[r[2] for r in example]
    beta, B, D, C=[sm[key] for key in ('BETA','B','D','C')]
    pi=arb.pi()
    coeff1,coeff2,coeffY=1/(2*pi),1/(8*pi**2),1/(32*pi**3)
    start=[*initial[:2],1/arb('.117').union(arb('.119'))]
    log0=[arb(v).log() for v in ('.98','.016','.010')]
    max_b=[beta[i]+16*(one[0][i]+one[1][i]) for i in range(3)]
    max_B=[[B[i][j]+16*(two[0][i][j]+two[1][i][j]) for j in range(3)] for i in range(3)]
    # 任意正规范解直到 L=29.4 都满足 y_f^2<2。三个逆耦合的上屏障
    # 60,30+.6L,9+1.2L 的外向导数为负；对汤川用合作的对数上比较系统。
    assert all(v>=0 for row in one for v in row) and all(v>=0 for m in two for row in m for v in row)
    assert all(v>=0 for i,row in enumerate(B) for j,v in enumerate(row) if (i,j)!=(2,2))
    assert all(0<u<v for u,v in zip(start,(60,30,9)))
    for i,slope in enumerate((0,arb('.6'),arb('1.2'))):
        assert -beta[i]*coeff1+max(0,-B[i][i])/(9)*coeff2+2*sum(D[i])*coeffY<slope
    upper=log0[:]
    for n in range(295):
        dt=arb('.1')
        time=(n+1)*dt
        inv=[1/arb(60),1/(30+arb('.6')*time),1/(9+arb('1.2')*time)]
        yy=[(2*(u+dt/3)).exp() for u in upper]
        assert max(yy)<2
        h=3*yy[0]+3*yy[1]+yy[2]
        rhs=[((arb('4.5')*yy[0]+arb('1.5')*yy[1]+yy[2])-4*pi*dot(C[0],inv))/(16*pi**2),
             ((arb('1.5')*yy[0]+arb('4.5')*yy[1]+yy[2])-4*pi*dot(C[1],inv))/(16*pi**2),
             ((h+arb('1.5')*yy[2])-4*pi*dot(C[2],inv))/(16*pi**2)]
        assert max(rhs)<arb(1)/3
        upper=[(u+dt*f).upper() for u,f in zip(upper,rhs)]
    assert max_b[1]-max_b[0]>27
    assert sum(max_B[1][:2])-sum(max_B[0][:2])>0 and max_B[1][2]>max_B[0][2]
    assert (max_b[1]-max_b[0])*coeff1-2*sum(abs(v-w) for v,w in zip(D[0],D[1]))*coeffY>0
    # 阈值全部开启后，在 X1=X2 的边界上 (X1-X2)' 严格为正。
    # 已进入 X1>X2 的轨迹可以提前终止，包括其正耦合解余下的整个区间。
    projections=((1,-1,0),(1,0,-1),(0,1,-1),(-1,2,-1),(2,-1,-1))
    contrasts=[(dot(v,start),dot(v,beta)*coeff1,[dot(v,r)*coeff1 for r in one],v) for v in projections]
    dt=arb(1)/128
    steps=3764
    assert steps*dt>limit
    knots=[times[q] for q in groups]
    def step(n,shift,state):
        corr,logs=state
        left,right=n*dt,(n+1)*dt
        t=left.union(right)
        if shift is None:
            opened=elapsed=[arb(0),arb(0)]
        else:
            low,high=left+shift.lower(),right+shift.upper()
            opened=[arb(sum(v<low for v in group)).union(sum(not(v>high) for v in group)) for group in knots]
            # 阈值早于 M_Z 时从 L=0 起累计；持续时间 min(t,(t+s-v)_+) 对 t,s 递增。
            elapsed=[sum(left.min((low-v).max(0)) for v in group).union(
                     sum(right.min((high-v).max(0)) for v in group)) for group in knots]
        x=[start[i]+corr[i]-(beta[i]*t+one[0][i]*elapsed[0]+one[1][i]*elapsed[1])*coeff1
           +arb(0,dt) for i in range(3)]
        if any(v<3 or v.contains(0) for v in x):
            return None,False,False
        inv=[1/v for v in x]
        yy=[(2*(v+arb(0,dt/3))).exp() for v in logs]
        coeff=[[B[i][j]+opened[0]*two[0][i][j]+opened[1]*two[1][i][j] for j in range(3)] for i in range(3)]
        dx=[-dot(row,inv)*coeff2+dot(r,yy)*coeffY for row,r in zip(coeff,D)]
        h=3*yy[0]+3*yy[1]+yy[2]
        dy=[(arb('1.5')*v+h-4*pi*dot(row,inv))/(16*pi**2) for row,v in zip(C,(yy[0]-yy[1],yy[1]-yy[0],yy[2]))]
        if any(not abs(v)<1 for v in dx) or any(not abs(v)<arb(1)/3 for v in dy):
            return None,False,False
        clear=False
        for idx,(base,b,m,vec) in enumerate(contrasts):
            residual=base+dot(corr,vec)-b*t-dot(m,elapsed)+arb(0).union(dt*dot(dx,vec))
            if not residual.contains(0):
                if idx==0 and shift is not None and left+shift.lower()>last and residual>0:
                    return None,True,True
                clear=True
                break
        return ([u+dt*v for u,v in zip(corr,dx)],[u+dt*v for u,v in zip(logs,dy)]),clear,False

    plain=[([arb(0)]*3,log0)]
    for n in range(steps):
        state,clear,_=step(n,None,plain[-1])
        assert state is not None and clear,('SM',n)
        plain.append(state)
    def exclude(lo,hi):
        sh=lo.union(hi)
        # 仅为节省积分工作；用严格不等式核准已有 SM 前缀可复用。
        index=max(0,min(steps,int(float((first-hi)/dt))-1))
        assert index==0 or index*dt+sh.upper()<first
        state=plain[index]
        for n in range(index,steps):
            state,clear,barrier=step(n,sh,state)
            if barrier:
                return True
            if not clear or state is None:
                return False
        return True
    left=2*(arb(701)/724).log()
    right=2*(arb(757)/724).log()
    # 所有阈值起于 M_Z 以上正是原尖阈值处方的入口条件。
    ranges=[(arb('-3.5'),left),(right,arb('27.013')),(arb('27.013'),arb(30))]
    assert start[0]>start[1] and last<30
    assert arb('-3.5')<2*(arb(128)/724).log()
    for lo,hi in ranges:
        n=max(1,int(float(hi-lo)*2)+1)
        pending=[(lo+(hi-lo)*i/n,lo+(hi-lo)*(i+1)/n) for i in range(n)]
        leaves=0
        while pending:
            l,r=pending.pop()
            if exclude(l,r):
                leaves+=1
            else:
                assert r-l>arb('.00002'),(l,r)
                mid=(l+r)/2
                pending.extend(((l,mid),(mid,r)))
        print('排除的阈值平移区间与末段数',lo.str(12),hi.str(12),leaves,flush=True)
    # 显式检查离合格分母最近的左右两个额外曲率整数对。
    for k,d in (((30,701),(35,757)) if neighbors else ()):
        effective=arb(k)*(arb(724)/d)**2
        curve=nominal_curve*(32/effective)**4*(1-arb(0).union(arb('1e-57')))
        vacuum=joint['vacuum']
        scaled={charge:[v*(arb(724)/d)**2*(arb(k)/32 if charge%2==0 else 1) for v in band]
                for charge,band in layers.items()}
        fields=vacuum['minimal']['fields'](scaled)
        residual=[(w+(arb(-g)/6 if g<0 else arb(g)/12),z) for g,w,z in fields]
        residual+=vacuum['low_sector'](4,a*(32/effective)**2)
        assert vacuum['newton'](residual,arb('.03'))<1<vacuum['newton'](residual,arb('.033'))
        assert arb('.998')<curve<arb('1.002')
        with redirect_stdout(StringIO()):
            energy,_,strong=gauge['crossing_enclosure'](model,initial,layers,False,d,'29.04','.1184',recenter=True)
        assert (strong>arb('.119') if d==701 else strong<arb('.117'))
        assert all(k*unit*(arb(724)/d)**2*v>energy for v in even)
        print('邻近整数、曲率和两圈交点强耦合',k,d,curve.str(12),strong.lower().str(10),strong.upper().str(10))
    print('全正整数只余',matching,'偶阈值对 L 的严格上界',limit.str(13))


if __name__=='__main__':
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)

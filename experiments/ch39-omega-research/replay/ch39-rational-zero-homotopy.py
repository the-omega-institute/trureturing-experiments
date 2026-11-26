"""有理圆点到两枚实测零点的荷负四正锥延拓及半径九处的四荷谱隙。"""
from pathlib import Path
from runpy import run_path
import sys
import numpy as np
from flint import acb, acb_mat, arb, arb_mat, ctx

CHARGES=(-4,1,2,-3,6)
I=acb_mat([[int(i==j) for j in range(19)] for i in range(19)])
I2=arb_mat([[1,0],[0,1]])


def adj(M):
    return M.conjugate().transpose()


def column(z):
    return arb_mat([[z.real],[z.imag]])


def make_models(root,local):
    old=run_path(str(root/'ch39-two-index-chord.py'))
    geometry=run_path(str(root/'ch39-unlabeled-native-connection.py'))
    rational=run_path(str(local))
    shell=run_path(str(root/'ch39-three-ray-chiral-boundary.py'))
    ctx.prec=192
    ends=rational['algebraic_ends'](geometry)
    models={}
    for q in CHARGES:
        A,D,B=ends(q)
        X,Y,C=old['ends'](shell,q,False)
        assert all((b-c).contains(0) for b,c in zip(B.entries(),C.entries()))
        models[q]=(A,D,X-A,Y-D,B)
    return models


def center_root(model,p,guess):
    A,D,P,S,B=model
    H=A+p*P
    d=D+p*S
    t,e=guess
    for _ in range(5):
        V=(H+t*d-e*I).solve(B)
        G=adj(B)*V
        tpart=-(adj(V)*d*V)[0,1]
        epart=(adj(V)*V)[0,1]
        J=arb_mat([[tpart.real,epart.real],[tpart.imag,epart.imag]])
        step=J.solve(column(G[0,1]))
        t,e=(t-step[0,0]).mid(),(e-step[1,0]).mid()
    return t,e


def root_cell(model,p,h,guess):
    A,D,P,S,B=model
    tc,ec=center_root(model,p,guess)
    Dc=D+p*S
    T,E=tc+arb(0,h*arb('.15')),ec+arb(0,h*arb('.06'))
    Pc=P+tc*S
    z=arb(0,h)
    M=A+p*P+tc*Dc-ec*I
    V0=M.solve(B)
    G0=adj(B)*V0
    tj=-(adj(V0)*Dc*V0)[0,1]
    ej=(adj(V0)*V0)[0,1]
    pre=arb_mat([[tj.real,ej.real],[tj.imag,ej.imag]]).inv()
    x,y=T-tc,E-ec
    N=M.inv()
    def fn(m):
        return sum(abs(v).upper()**2 for v in m.entries()).sqrt().upper()
    small=(x,z,x*z,-y)
    mats=(N*Dc,N*Pc,N*S,N)
    eta=sum(abs(v).upper()*fn(C) for v,C in zip(small,mats))
    if not eta<arb('.6'):
        raise ArithmeticError("预乘逆")
    first=acb_mat(19,2)
    for v,C in zip(small,mats):
        first+=v*(C*V0)
    rem=fn(first)*eta/(1-eta)
    unit=acb(arb(0,rem.upper()),arb(0,rem.upper()))
    V=V0-first+acb_mat([[unit]*2 for _ in range(19)])
    def parts(V):
        star=adj(V)
        return -(star*(Dc+z*S)*V),star*V,-(star*(Pc+x*S)*V)
    Tj,Ej,Pj=parts(V)
    J=arb_mat([[Tj[0,1].real,Ej[0,1].real],[Tj[0,1].imag,Ej[0,1].imag]])
    jdiff=I2-pre*J
    spread=(pre*column(Pj[0,1]))*z+(pre*column(G0[0,1]))
    # 映射区间与中心化根盒的算子范数均按两个轴的不同半宽计算。
    pad=(arb('.15')*h,arb('.06')*h)
    im=[];con=[]
    for i in range(2):
        err=sum(abs(jdiff[i,j]).upper()*pad[j] for j in range(2))
        con.append(err/pad[i])
        im.append(abs(spread[i,0]).upper()+err)
    if not all(v<s for v,s in zip(im,pad)) or not all(v<1 for v in con):
        raise ArithmeticError("根盒")
    # 一次均值求值同时利用三个 Hermitian 方向，避免对两个边权直接作宽区间求逆。
    tt,ee=tc+arb(0,im[0].upper()),ec+arb(0,im[1].upper())
    r=[]
    for j in range(2):
        g=G0[j,j].real+Tj[j,j].real*(tt-tc)+Ej[j,j].real*(ee-ec)+Pj[j,j].real*z
        r.append(-1/g)
    if not -1<r[0]<0<r[1]:
        raise ArithmeticError("边权")
    return (tc,ec), (T,E), (tt,ee), r, max(c.upper() for c in con)


def spectator(model,p,h,t,w):
    A,D,P,S,B=model
    tc=t.mid()
    wc=[a.mid() for a in w]
    midpoint=A+p*P+tc*(D+p*S)+B*acb_mat([[wc[0],0],[0,wc[1]]])*adj(B)-9*I
    center=np.array(midpoint.tolist(),complex)
    vals,vec=np.linalg.eigh(center)
    assert sum(vals<0)==3
    # 真矩阵在此数值合同基中用 Arb 重新计算，故基向量只作有理近似建议。
    V=acb_mat((vec/np.sqrt(abs(vals))).tolist())
    A0=adj(V)*midpoint*V
    for i in range(19):A0[i,i]-=1 if vals[i]>0 else -1
    td=t-tc
    P1=adj(V)*(P+tc*S)*V
    T1=adj(V)*(D+p*S)*V
    T2=adj(V)*S*V
    col=adj(V)*B
    err=A0+arb(0,h)*P1+td*T1+arb(0,h)*td*T2
    for j in range(2):
        v=acb_mat([[col[i,j]] for i in range(19)])
        err+=(w[j]-wc[j])*v*adj(v)
    norm=sum(abs(z).upper()**2 for z in err.entries()).sqrt()
    if not norm<1:
        raise ArithmeticError("半径九")
    # 谱到9距离至少 (1-||err||) / ||V||²，Frobenius 上界即可。
    gap=(1-norm)/sum(abs(z).upper()**2 for z in V.entries())
    return gap.lower()


def verify(root,local):
    models=make_models(root,local)
    def guess(p):
        return arb('.0198')+arb('.056')*p,arb('9.009')-arb('.01')*p
    def assess(p,h):
        centered,outer,inner,w,rho=root_cell(models[-4],p,h,guess(p))
        gaps=[spectator(models[q],p,h,inner[0],w) for q in CHARGES[1:]]
        if not min(gaps)>arb(".0004"):
            raise ArithmeticError("旁观谱隙下界")
        return centered,outer,inner,w,rho,min(gaps)
    records=[]
    def cover(a,b):
        p,h=(a+b)/2,(b-a)/2
        try: info=assess(p,h)
        except (ArithmeticError,ZeroDivisionError):
            if h<arb(2)**-18:raise
            cover(a,p)
            cover(p,b)
            return
        records.append((a,b,*info))
    cover(arb(0),arb(1))
    # 在公共端点单独圈定根，验证两个单元的唯一根一致。
    for k in range(len(records)+1):
        p=records[k][0] if k<len(records) else arb(1)
        _,_,inner,_,_=root_cell(models[-4],p,arb('1e-12'),guess(p))
        if k in (0,len(records)):
            anchor=(".019827825195380612","9.0078347438573692") if k==0 else (".070974807554382935","8.9989005066886309")
            assert all((arb(ref)+arb(0,arb("1e-12"))).contains(z) for ref,z in zip(anchor,inner))
        for j in (k-1,k):
            if 0<=j<len(records):
                assert all(z.contains(v) for z,v in zip(records[j][3],inner))
    print("根盒数量、最大预乘压缩率、半径九的统一隙下界",len(records),
          max(rec[6] for rec in records).str(12),min(rec[7] for rec in records).str(12))
    for i,label in enumerate(("弦坐标","双根能量","第一边电导","第二边电导")):
        values=[rec[4][i] if i<2 else 1+rec[5][i-2] for rec in records]
        print(label,min(z.lower() for z in values).str(12),max(z.upper() for z in values).str(12))


if __name__=='__main__':
    here=Path(__file__).resolve().parent
    root=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else here
    del sys.argv[1:]
    verify(root,root/'ch39-algebraic-two-clocks.py')

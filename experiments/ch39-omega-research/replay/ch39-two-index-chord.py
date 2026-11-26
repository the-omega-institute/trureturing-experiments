"""第36、1枚实际零点的壳矩阵弦在带颜色重数的荷 −4 部门形成正向 Weyl 点。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, arb_mat, ctx


EDGES = ((0,2),(5,14))
CHARGES = (1,-4,2,-3,6)


def adj(m):
    return m.conjugate().transpose()


def frobenius(m):
    return sum(abs(z).upper()**2 for z in m.entries()).sqrt()


def ends(shell, q, flat):
    result=[]
    for n in (36,1):
        h=acb.zeta_zero(n).imag
        H=shell["laplacian"](q,h)
        if flat:
            soft=shell["soft"]
            delta=arb(23)/25*arb(2)**-20/(h*h+soft["point"]["study"]["offset"])
            for j in soft["star"]["orbit"]:
                i=soft["pole"]
                H[i,i]+=delta
                H[j,j]+=delta
                H[i,j]-=delta
                H[j,i]-=delta
        result.append(H)
    A,Z=result
    B=acb_mat(19,2)
    for k,(i,j) in enumerate(EDGES):
        assert abs(A[i,j]).contains(1) and A[i,j].overlaps(Z[i,j])
        B[i,k]=1
        B[j,k]=A[j,i]
    return A,Z-A,B


def green(A,D,B,eye,t,e):
    V=(A+t*D-e*eye).solve(B)
    G=adj(B)*V
    x=-(adj(V)*D*V)[0,1]
    y=(adj(V)*V)[0,1]
    J=arb_mat([[x.real,y.real],[x.imag,y.imag]])
    return G,J,V


def gap_and_count(matrix,energy):
    vals=matrix.eig(multiple=True)
    assert len(vals)==19 and all(z.imag.contains(0) for z in vals)
    assert all(not z.real.overlaps(energy) for z in vals)
    return sum(z.real<energy for z in vals),min(abs(z.real-energy).lower() for z in vals)


def tetrahedral_degree(A, D, B, t, inc):
    # 此顶点序对 (t,w1,w2) 为负向，故取标准四面体边界的反向。
    directions=((1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1))
    faces=((1,3,2),(0,2,3),(0,3,1),(0,1,2))
    scale=arb("1e-13")
    frames=[]
    for h,x,y in directions:
        weights=acb_mat([[inc[0].mid()+scale*x,0],[0,inc[1].mid()+scale*y]])
        M=A+(t.mid()+scale*h)*D+B*weights*adj(B)
        vals,vec=M.eig(right=True,algorithm="rump")
        # 任由相近根在双精度下并列，使用区间实部比较。
        order=sorted(range(19),key=lambda j:vals[j].real)
        assert all(vals[i].real<vals[j].real for i,j in zip(order[:3],order[1:4]))
        cols=[]
        for j in order[:3]:
            z=[vec[k,j] for k in range(19)]
            n=sum(abs(a)**2 for a in z).sqrt()
            cols.append([a/n for a in z])
        frames.append(acb_mat(cols).transpose())
    angles=[]
    for i,j,k in faces:
        a=(adj(frames[j])*frames[i]).det()
        b=(adj(frames[k])*frames[j]).det()
        c=(adj(frames[i])*frames[k]).det()
        assert abs(a)>0 and abs(b)>0 and abs(c)>0
        angle=(a*b*c).arg()
        assert abs(angle)<arb.pi()
        angles.append(angle)
    degree=(sum(angles)/(2*arb.pi())).unique_fmpz()
    assert degree is not None
    return int(degree),angles


def certify(shell,flat):
    A,D,B=ends(shell,-4,flat)
    eye=acb_mat([[int(i==j) for j in range(19)] for i in range(19)])
    eye2=arb_mat([[1,0],[0,1]])
    t0,e0=arb(".070974807555"),arb("8.99890050669")
    for _ in range(2):
        G,J,_=green(A,D,B,eye,t0,e0)
        y=J.solve(arb_mat([[G[0,1].real],[G[0,1].imag]]))
        t0,e0=(t0-y[0,0]).mid(),(e0-y[1,0]).mid()
    G0,J0,_=green(A,D,B,eye,t0,e0)
    pre=J0.inv()
    residual=max(abs(z).upper() for z in (pre*arb_mat([[G0[0,1].real],[G0[0,1].imag]])).entries())
    for radius in map(arb,("1e-10","1e-23")):
        t,e=t0+arb(0,radius),e0+arb(0,radius)
        G,J,V=green(A,D,B,eye,t,e)
        defect=eye2-pre*J
        rho=max(sum(abs(defect[i,j]) for j in range(2)).upper() for i in range(2))
        assert rho<1 and residual+rho*radius<radius
    inc=[-1/G[j,j].real for j in range(2)]
    assert -1<inc[0]<0<inc[1] and 0<t<1
    below,gap=gap_and_count(A+t0*D,e0)
    assert below==3 and gap>(1+frobenius(D))*radius
    v=acb_mat([[V[i,0]] for i in range(19)])
    v=v/(adj(v)*v)[0,0].real.sqrt()
    w=acb_mat([[V[i,1]] for i in range(19)])
    w=w-v*(adj(v)*w)[0,0]
    w=w/(adj(w)*w)[0,0].real.sqrt()
    F=acb_mat([[v[i,0],w[i,0]] for i in range(19)])
    pauli=[]
    deriv=[D]
    for j in range(2):
        b=acb_mat([[B[i,j]] for i in range(19)])
        deriv.append(b*adj(b))
    for m in deriv:
        X=adj(F)*m*F
        pauli.append([X[0,1].real,-X[0,1].imag,(X[0,0].real-X[1,1].real)/2])
    speed=arb_mat(pauli).det()
    assert speed>arb("7e-6")
    print("端点软辐边", "单位权" if flat else "实际亏额")
    print("t、双根能量、两条正电导", t.str(18),e.str(18),*[(1+z).str(18) for z in inc])
    print("根盒预乘压缩率上界、Pauli 速度行列式",rho.str(9),speed.str(14))
    degrees=[]
    spectators=[]
    third_levels, fourth_levels = [], []
    for q in CHARGES:
        M,C,Q=ends(shell,q,flat)
        assert frobenius(C)<20
        degree,angles=tetrahedral_degree(M,C,Q,t,inc)
        degrees.append(degree)
        if q==-4:
            print("外向小四面体四面相角",*[z.str(12) for z in angles])
            continue
        T=arb_mat([[inc[0],0],[0,inc[1]]])
        L=M+t*C+Q*acb_mat(T)*adj(Q)
        center=acb_mat([[L[i,j].mid() for j in range(19)] for i in range(19)])
        center=(center+adj(center))/2
        count,gap=gap_and_count(center,e0)
        certified=gap-frobenius(L-center)-radius
        assert count==3 and certified>arb(".009")
        ordered = sorted((z.real for z in center.eig(multiple=True)), key=float)
        perturbation = arb(0, frobenius(L-center).upper())
        third_levels.append(ordered[2]+perturbation)
        fourth_levels.append(ordered[3]+perturbation)
        spectators.append((q,count))
        print("旁观荷、交点能量以下态数、距该能量严格下界、722 态负谱丛第二陈数",
              q,count,certified.lower().str(12),count)
    assert spectators==[(1,3),(2,3),(-3,3),(6,3)]
    low = max(third_levels, key=float)
    high = min(fourth_levels, key=float)
    assert all(v < 9 for v in third_levels) and all(9 < v for v in fourth_levels)
    assert low > 8 and high < 10
    assert all(low > v for v in third_levels if v is not low)
    assert all(high < v for v in fourth_levels if v is not high)
    print("四种配对共用的半径开区间端点与唯一整数半径", low.str(14), high.str(14), 9)
    index=sum(n*d*d for n,d in zip((6,3,3,2,1),degrees))
    assert degrees==[0,1,0,0,0] and index==3
    print("五荷球面第一陈数与标准重数的同荷外积指数",degrees,index)


def verify(root):
    shell=run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec=256
    certify(shell,True)
    certify(shell,False)


if __name__=="__main__":
    root=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
    del sys.argv[1:]
    verify(root)

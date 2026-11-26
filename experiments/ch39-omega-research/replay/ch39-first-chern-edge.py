"""原生双钟连续高度与两条正边权形成的第三、第四带 Weyl 点。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, arb_mat, ctx


EDGES = ((0,2),(5,14))
Q = -3


def adj(m):
    return m.conjugate().transpose()


def calculate(shell, height, energy):
    H = shell["laplacian"](Q,height)
    dH = acb_mat(19,19)
    B, dB = acb_mat(19,2), acb_mat(19,2)
    geo, soft = shell["geo"], shell["soft"]
    offset = soft["point"]["study"]["offset"]
    delta_prime = -2*height * arb(23)/25 * arb(2)**-20 / (height**2+offset)**2
    for i,j,d,_ in geo["edges"]:
        softened = soft["pole"] in (i,j) and (i in soft["star"]["orbit"] or j in soft["star"]["orbit"])
        wp = -delta_prime if softened else arb(0)
        w = abs(H[i,j])
        freq = geo["theta"](d)-geo["theta"](244) if d in (724,804) else arb(0)
        dH[i,i] += wp
        dH[j,j] += wp
        dH[i,j] = (wp/w+acb(0,Q*freq))*H[i,j]
        dH[j,i] = dH[i,j].conjugate()
        if (i,j) in EDGES:
            k = EDGES.index((i,j))
            assert not softened and d not in (724,804) and w.contains(1)
            B[i,k] = 1
            B[j,k] = H[i,j].conjugate()
            dB[j,k] = dH[i,j].conjugate()
    eye = acb_mat([[int(i==j) for j in range(19)] for i in range(19)])
    V = (H-energy*eye).solve(B)
    G = adj(B)*V
    T = adj(dB)*V
    ax = (T+adj(T)-adj(V)*dH*V)[0,1]
    ay = (adj(V)*V)[0,1]
    jac = arb_mat([[ax.real,ay.real],[ax.imag,ay.imag]])
    return G, jac, (H,dH,B,dB,V)


def verify(root):
    shell=run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec=256
    h0=acb.zeta_zero(22).imag
    s0=arb("-1.98195516305233")
    e0=arb("8.92931408205483")
    original_s, original_e = s0,e0
    for radius in map(arb, ("1e-11","1e-25")):
        G0,J0,_=calculate(shell,h0+s0,e0)
        R=J0.inv()
        y=R*arb_mat([[G0[0,1].real],[G0[0,1].imag]])
        s=s0+arb(0,radius)
        e=e0+arb(0,radius)
        G,J,(H,D,B,dB,V)=calculate(shell,h0+s,e)
        eye=arb_mat([[1,0],[0,1]])
        err=eye-R*J
        contraction=max(sum(abs(err[i,j]) for j in range(2)).upper() for i in range(2))
        residual=max(abs(v).upper() for v in y.entries())
        assert contraction<1 and residual+radius*contraction<radius
        if radius>arb("1e-20"):
            s0=(s0-y[0,0]).mid()
            e0=(e0-y[1,0]).mid()
    assert abs(s-original_s)<arb("1e-12") and abs(e-original_e)<arb("1e-12")
    assert acb.zeta_zero(21).imag<h0+s<h0
    # 钟边的导数算子范数逐条加总，另有三条软辐边。
    g=shell["geo"]
    freq=[abs(Q*(g["theta"](d)-g["theta"](244))) for _,_,d,_ in g["edges"] if d in (724,804)]
    assert len(freq)==12 and sum(freq)+1<20
    assert h0+s>80 and shell["soft"]["point"]["study"]["offset"]>0
    t=[-1/G[j,j].real for j in range(2)]
    assert -1<t[0]<0<t[1]
    H0=calculate(shell,h0+s0,e0)[2][0]
    ev=H0.eig(multiple=True)
    # 用 Frobenius 范数控制整个小盒中的本征值移动。
    change=sum(abs(z).upper()**2 for z in (H-H0).entries()).sqrt()
    assert all(abs(z.real-e0)>change+radius for z in ev) and sum(z.real<e0 for z in ev)==3
    # 在精确重根的二维核中用固定顺序的 Gram–Schmidt 正交化。
    v=acb_mat([[V[i,0]] for i in range(19)])
    v=v/(adj(v)*v)[0,0].real.sqrt()
    w=acb_mat([[V[i,1]] for i in range(19)])
    w=w-v*(adj(v)*w)[0,0]
    w=w/(adj(w)*w)[0,0].real.sqrt()
    F=acb_mat([[v[i,0],w[i,0]] for i in range(19)])
    ders=[]
    for k in range(2):
        b=acb_mat([[B[i,k]] for i in range(19)])
        db=acb_mat([[dB[i,k]] for i in range(19)])
        D+=t[k]*(db*adj(b)+b*adj(db))
        ders.append(b*adj(b))
    pauli=[]
    for m in (D,*ders):
        sM=adj(F)*m*F
        pauli.append([sM[0,1].real,-sM[0,1].imag,(sM[0,0].real-sM[1,1].real)/2])
    det=arb_mat(pauli).det()
    assert det > arb("3.45e-5")
    print("高度 h、第22枚相对位移、交叉本征值", (h0+s).str(20),s.str(20),e.str(20))
    print("两原边及其正电导",EDGES,*[(1+z).str(19) for z in t])
    print("零方程的雅可比行列式、压缩率上界、Pauli 速度行列式",J.det().str(10),contraction.str(10),det.str(12))
    print("按 (h,w1,w2) 定向，采用 i Tr(P dP∧dP)/(2π) 时第三带下支的陈数",1 if det>0 else -1)


if __name__=="__main__":
    root=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
    del sys.argv[1:]
    verify(root)

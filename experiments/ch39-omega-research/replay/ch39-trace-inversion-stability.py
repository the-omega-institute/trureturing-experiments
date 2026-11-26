"""荷一两轨道三阶迹反推双钟的局部误差放大。"""
from math import factorial
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, arb_mat, ctx


NUMBER = 2_306_888_100_085


def column(v):
    return arb_mat([[x] for x in v])


def verify(root):
    exact = run_path(str(root/"ch39-mincut-real-triangle-tomography.py"))
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 448
    h = [z.imag for z in acb.zeta_zeros(NUMBER,8)]
    x = [arb(0),*[y-h[0] for y in h[1:]]]
    assert len(x)==8 and all(l<r for l,r in zip(x,x[1:])) and x[-1]<2
    print("八个零点距首点的高度差",*[v.str(16) for v in x])
    angles = [(1-arb(d)/1734).acos() for d in (244,724,804)]
    w = [a-angles[0] for a in angles[1:]]
    amplitudes = [-36*exact["unit"](label).imag for label in ((510,1),(1173,-2))]
    A = [amp*(freq*h[0]).sin() for amp,freq in zip(amplitudes,w)]
    B = [amp*(freq*h[0]).cos() for amp,freq in zip(amplitudes,w)]
    print("两角频率、中心余弦和正弦系数",*[v.str(14) for v in (*w,*A,*B)])
    def jets(freq,k):
        return [[t**k*((freq*t).cos() if k%2==0 else -(freq*t).sin())*(-1)**(k//2) for t in x],
                [t**k*((freq*t).sin() if k%2==0 else (freq*t).cos())*(-1)**(k//2) for t in x]]
    J = []
    for j,freq in enumerate(w):
        c,s = jets(freq,0)
        dc,ds = jets(freq,1)
        J.extend([c,s,[A[j]*a+B[j]*b for a,b in zip(dc,ds)]])
    J = arb_mat(J).transpose()
    gram = J.transpose()*J
    eigen = sorted(gram.eig(),key=lambda z:float(z.real))
    assert all(z.imag.contains(0) and z.real>0 for z in eigen)
    assert len(eigen)==6 and all(a.real<b.real for a,b in zip(eigen,eigen[1:]))
    singular = [z.real.sqrt() for z in eigen]
    print("奇异值由小到大",*[v.str(16) for v in singular])
    P = gram.solve(J.transpose())
    print("各参数欧氏噪声放大",*[(sum(P[i,j]**2 for j in range(8))).sqrt().str(15) for i in range(6)])
    rowsums = [sum(abs(P[i,j]) for j in range(8)) for i in range(6)]
    print("无穷范数各行放大",*[v.str(15) for v in rowsums])
    print("条件数",(singular[-1]/singular[0]).str(15))
    assert arb("3.40294066e-8")<singular[0]<arb("3.40294067e-8")
    assert arb("6.97217860")<singular[-1]<arb("6.97217861")
    print("欧氏最坏反演放大",(1/singular[0]).str(16))
    order = 6
    tensor = []
    for freq in w:
        tensor.append([[P*column(values) for values in jets(freq,k)] for k in range(order+2)])

    def derivative_bound(radius):
        # 对已乘 P 的列作泰勒展开，避免在邻域上丢失中心处的抵消。
        total = [arb(0) for _ in range(6)]
        for j in range(2):
            for i in range(6):
                def variation(n,ca,cb):
                    changes = sum(abs(ca*tensor[j][n+k][0][i,0]+cb*tensor[j][n+k][1][i,0])*
                                  radius**k/factorial(k) for k in range(1,order+1))
                    tail = (abs(ca)+abs(cb))*rowsums[i]*x[-1]**(n+order+1)*radius**(order+1)/factorial(order+1)
                    return changes+tail
                freq_shift = variation(1,A[j],B[j])
                coeff_shift = sum(abs(tensor[j][1][k][i,0])+variation(1,1-k,k) for k in range(2))
                total[i] += variation(0,1,0)+variation(0,0,1)+freq_shift+radius*coeff_shift
        return max(v.upper() for v in total)

    radius = arb(1)/8000
    change = derivative_bound(radius)
    pnorm = max(v.upper() for v in rowsums)
    amplification = pnorm/(1-change)
    noise = radius/amplification
    print("非线性参数立方半径、预乘后偏差、误差放大和可容许扰动",radius.str(10),
          change.str(16),amplification.str(16),noise.str(16))
    assert change<arb(".434") and amplification<arb("8.67e7") and noise>arb("1.44e-12")
    p0 = [A[0],B[0],w[0],A[1],B[1],w[1]]
    def signal(p):
        return column([sum(p[k]*(p[k+2]*t).cos()+p[k+1]*(p[k+2]*t).sin() for k in (0,3)) for t in x])
    y0 = signal(p0)
    assert arb(".27215586483699")<y0[0,0]<arb(".27215586483700")
    for j,raw in enumerate(h):
        assert y0[j,0].overlaps(sum(a*(f*raw).sin() for a,f in zip(amplitudes,w)))
    # 取第二钟正弦系数对噪声最敏感的方向，迭代仅用来生成一个可直接核验的扰动。
    k = 4
    norm = sum(P[k,i]**2 for i in range(8)).sqrt()
    target = y0+column([P[k,i]/norm*arb("5e-12") for i in range(8)])
    alt = p0[:]
    for _ in range(7):
        dp = P*(target-signal(alt))
        alt = [(v+dp[i,0]).mid() for i,v in enumerate(alt)]
    diff = signal(alt)-y0
    print("有限扰动的六个参数变化",*[(v-u).str(14) for u,v in zip(p0,alt)])
    largest = max(abs(z).upper() for z in diff.entries())
    assert max(abs(v-u).upper() for u,v in zip(p0,alt))<radius
    assert arb("1.8853e-6")<alt[2]-p0[2]<arb("1.8854e-6")
    assert arb("-1.2220e-6")<alt[5]-p0[5]<arb("-1.2219e-6")
    assert max(abs(z).lower() for z in diff.entries())>arb("3.3049e-12")
    assert largest<arb("3.305e-12")
    offset = shell["soft"]["point"]["study"]["offset"]
    raw = [z*arb(23)/25*arb(2)**-20/(v*v+offset) for z,v in zip(diff.entries(),h)]
    raw_error = max(abs(z).upper() for z in raw)
    assert raw_error<arb("8.14e-42")
    print("八个归一化读数的最大与欧氏变化、未归一化迹差的最大变化",
          largest.str(15),sum(z*z for z in diff.entries()).sqrt().str(15),raw_error.str(15))



if __name__=="__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)

"""以四枚实际正零点补充原八点，验证六参数反演的稳定性。"""
from math import factorial

from flint import acb, arb, arb_mat, ctx


N = 2_306_888_100_085
EXTRA = (-279,-232,-45,99)


def rows(times,p):
    columns = []
    for a,b,w in (p[:3],p[3:]):
        trig = [(w*t).sin_cos() for t in times]
        columns.extend([[c for s,c in trig],[s for s,c in trig],
                        [t*(-a*s+b*c) for t,(s,c) in zip(times,trig)]])
    return arb_mat(columns).transpose()


def spectrum(jacobian):
    gram = jacobian.transpose()*jacobian
    eig = sorted(gram.eig(),key=lambda z:float(z.real))
    assert len(eig)==6 and all(z.imag.contains(0) and z.real>0 for z in eig)
    assert all(a.real<b.real for a,b in zip(eig,eig[1:]))
    return [z.real.sqrt() for z in eig],gram.solve(jacobian.transpose())


def verify():
    ctx.prec=256
    heights=[z.imag for z in acb.zeta_zeros(N,8)]
    base=heights[0]
    times=[arb(0),*[y-base for y in heights[1:]]]
    assert all(a<b for a,b in zip(times,times[1:])) and times[-1]<2
    more=[acb.zeta_zero(N+k).imag-base for k in EXTRA]
    print("四个附加零点的编号偏移、高度差",EXTRA,*[z.str(16) for z in more],flush=True)
    theta=[(1-arb(d)/1734).acos() for d in (244,724,804)]
    p=[]
    for angle,(r,i) in zip(theta[1:],((510,1),(1173,-2))):
        w=angle-theta[0]
        z=acb(r,i*arb(867).sqrt())
        coefficient=-36*(z/abs(z)).imag
        s,c=(w*base).sin_cos()
        p.extend((coefficient*s,coefficient*c,w))
    short,_=spectrum(rows(times,p))
    times+=more
    J=rows(times,p)
    sig,P=spectrum(J)
    print("原八点奇异值",*[s.str(13) for s in short])
    print("扩充后奇异值",*[s.str(13) for s in sig])
    assert sig[0]>arb(".8573") and short[3]<arb(".171")
    assert sig[0]/short[0]>25_190_000
    print("线性改善倍数、十二点条件数",(sig[0]/short[0]).str(13),(sig[-1]/sig[0]).str(13))
    print("频率行欧氏范数、绝对行和",*[(sum(P[k,j]**2 for j in range(12)).sqrt().str(13),
                                   sum(abs(P[k,j]) for j in range(12)).str(13)) for k in (2,5)])
    print("无穷算子范数",max(sum(abs(P[k,j]) for j in range(12)).upper() for k in range(6)).str(13))
    # 在同一个小立方体内，先左乘 P，再按每项频率的偏移展开导数。
    maxt=max(abs(t).upper() for t in times)
    degree=6
    jets=[]
    for k in (0,3):
        band=[]
        for n in range(degree+2):
            phase=n*arb.pi()/2
            values=[(p[k+2]*t+phase).sin_cos() for t in times]
            band.append([P*arb_mat([[t**n*pair[c]] for t,pair in zip(times,values)]) for c in (1,0)])
        jets.append(band)
    def bounds(radius):
        D=arb_mat(6,6)
        for block,k in enumerate((0,3)):
            band=jets[block]
            for i in range(6):
                norm=sum(abs(P[i,j]) for j in range(12))
                def variation(n,a,b):
                    main=sum(abs(a*band[n+j][0][i,0]+b*band[n+j][1][i,0])*radius**j/factorial(j)
                             for j in range(1,degree+1))
                    tail=(abs(a)+abs(b))*norm*maxt**(n+degree+1)*radius**(degree+1)/factorial(degree+1)
                    return (main+tail).upper()
                D[i,k]=variation(0,1,0)
                D[i,k+1]=variation(0,0,1)
                D[i,k+2]=variation(1,p[k],p[k+1])+radius*sum(
                    abs(band[1][c][i,0])+variation(1,1-c,c) for c in (0,1))
        return arb_mat([[D[i,j].upper() for j in range(6)] for i in range(6)])

    radius=arb(1)/8000
    d=bounds(radius)
    dnorm=max(sum(d[i,j] for j in range(6)).upper() for i in range(6))
    assert dnorm<arb(".856")
    eye=arb_mat([[int(i==j) for j in range(6)] for i in range(6)])
    row1=arb_mat([[sum(abs(P[i,j]) for j in range(12))] for i in range(6)])
    row2=arb_mat([[sum(P[i,j]**2 for j in range(12)).sqrt()] for i in range(6)])
    v=(eye-d).solve(row1)
    # |Δp| ≤ |P ΔF| + D |Δp|。对尚未估计的修正用逐分量的 Neumann 级数。
    correction=(eye-d).solve(d*row2)
    reverse=1/sig[0]+sum(z*z for z in correction.entries()).sqrt()
    safe=max(z.upper() for z in v.entries())
    print("整盒逐点噪声至六参数的放大上界",*[z.str(13) for z in v.entries()])
    print("欧氏逆放大上界、整盒最小奇异值下界、可容许逐点噪声",
          reverse.str(13),(1/reverse).str(13),(radius/safe).str(13))
    # ∞ 范数足够使重建迭代压缩。选 v_i ε 为每坐标半径则整盒映入自身。
    print("整盒预乘导数偏差的无穷范数上界",dnorm.str(13))
    assert reverse<arb("1.293") and 1/reverse>arb(".7738") and safe<arb("2.645")
    assert v[2,0]<arb(".2031") and v[5,0]<arb(".3637")
    assert radius/safe>arb("4.727e-5")
    # 有限扰动按十进制有理数指定；短、长基线直接重算同一对参数的观测。
    perturbations = [arb(v) for v in ("-5.1777326559295e-6","-0.00010340015819093","1.8853117842186e-6",
                                      "5.1777303306506e-6","0.00010412611599943","-1.2219050481320e-6")]
    assert all(abs(v)<radius for v in perturbations)
    def signal(t,param):
        return sum(param[k]*(t*param[k+2]).cos()+param[k+1]*(t*param[k+2]).sin() for k in (0,3))
    alternative = [a+b for a,b in zip(p,perturbations)]
    data = [signal(t,p) for t in times]
    differences = [signal(t,alternative)-y for t,y in zip(times,data)]
    assert arb(".27215586483699")<data[0]<arb(".27215586483700")
    print("四个新观测的原三阶迹差除以亏额",*[v.str(15) for v in data[8:]])
    print("同一有限扰动对八个原读数的最大绝对变化",max(abs(v).upper() for v in differences[:8]).str(14))
    print("对四个新读数的有符号变化",*[z.str(13) for z in differences[8:]])
    assert max(abs(z).upper() for z in differences[:8])<arb("3.306e-12")
    assert arb("-.000275051") < differences[8] < arb("-.000275050")


if __name__=="__main__":
    verify()

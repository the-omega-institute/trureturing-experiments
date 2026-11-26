"""中性参考固定时，荷依赖逐边电导对两座十六重层塔的必要反差与总和构造。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx


NUMBER=2_306_888_100_085


def verify(root):
    shell=run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec=256
    s0=arb(57528)/144905
    h=acb.zeta_zero(NUMBER).imag
    delta=arb(23)/25*arb(2)**-20/(h*h+shell["soft"]["point"]["study"]["offset"])
    c=15-delta/3
    spectrum={}
    for q in (1,2):
        roots=shell["laplacian"](q,h).eig(multiple=True)
        assert len(roots)==19 and all(z.imag.contains(0) and z.real>0 for z in roots)
        spectrum[q]=[z.real for z in roots]
    def heavy(q,t):
        folded=[t*z*abs(t*z-c) for z in spectrum[q]]
        low=sorted(z.lower() for z in folded)[3:]
        high=sorted(z.upper() for z in folded)[3:]
        return [a.union(b) for a,b in zip(low,high)]
    odd=heavy(1,arb(1))
    nominal=heavy(2,arb(1))
    target=32*sum(nominal)/sum(odd)
    def bound(t,s,neutral_drift):
        # L₂ 与同相位、荷一电导的矩阵之 Frobenius 距离受逐边反差控制。
        R=(s+18).sqrt()
        upper=t*R*((s+1).sqrt()+(t-1)*R+arb(19).sqrt()*neutral_drift)
        lower=s-18*neutral_drift
        assert lower>0
        return 19*upper/(16*lower)
    def bracket(f,a,b,n=48):
        assert f(a)<0<f(b)
        for _ in range(n):
            mid=(a+b)/2
            z=f(mid)
            if z<0:
                a=mid
            else:
                assert z>0
                b=mid
        return a,b
    native_lap=shell["laplacian"](0,h)
    sq=native_lap*native_lap
    s_native=sum(sq[i,i].real for i in range(19))/(c*c)-18
    cases=(("荷一即中性，任意参考",s0,False,arb("1.19148870")),
           ("三组逐边直径均受限",s0,True,arb("1.01019934")),
           ("荷一即中性，原参考",s_native,False,arb("1.58243783")))
    for label,s,drifting,floor in cases:
        def excess(t):
            return bound(t,s,t-1 if drifting else arb(0))-target
        top=1+s/20 if drifting else arb(2)
        a,b=bracket(excess,arb(1),top)
        assert a>floor
        print(label,"必要反差 t 的临界值及 log",a.union(b).str(16),a.union(b).log().str(16))
    def difference(t):
        return sum(heavy(2,t))-32*sum(nominal)
    a,b=bracket(difference,arb(2),arb(4))
    root=a.union(b)
    changed=heavy(2,root)
    eigen=sorted(spectrum[2],key=float)
    assert all(x<y for x,y in zip(eigen,eigen[1:]))
    assert a*eigen[3]>c and a*eigen[3]*(a*eigen[3]-c)>c*c/4
    assert arb("2.978020105")<root<arb("2.978020107")
    required=[]
    for lam,desired in zip(eigen[3:],nominal):
        assert 32*desired>c*c/4
        required.append((c+(c*c+128*desired).sqrt())/(2*lam))
    index=max(range(16),key=lambda j:float(required[j]))
    assert all(j==index or z<required[index] for j,z in enumerate(required))
    print("逐层重现的最强单层必要界来自第",index+4,"层，t、log(t)",
          required[index].str(16),required[index].log().str(16))
    print("固定原高度、原边相位的荷二同比放大根与 log",a.union(b).str(16),a.union(b).log().str(16))
    print("两塔目标总和比",target.str(16))
    for j in (0,15):
        print("第",j+4,"层实际、逐层目标及比例",changed[j].str(14),(32*nominal[j]).str(14),
              (changed[j]/(32*nominal[j])).str(14))
    # 中性参考与两荷均未锚定时，正电导相同、相位可独立的反例。
    edges={(i,j) for i,j,*_ in shell["geo"]["edges"]}
    triangle=next((i,j,k) for i,j in edges for k in range(j+1,19) if (i,k) in edges and (j,k) in edges)
    error=arb(2)*(len(edges)-3)/10**6
    ratio_floor=4/(16*(3+error)*error)
    assert len(edges)==135 and ratio_floor>target
    print("固定中性、两荷共享另一组权时的三角顶点和负通量比下界",triangle,ratio_floor.str(16))


if __name__=="__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)

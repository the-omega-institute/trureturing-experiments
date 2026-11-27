"""允许微扰顶 Yukawa 全区间后的独立中性规范阈值约束。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


ROOT = Path(__file__).resolve().parent


def verify():
    ctx.prec = 224
    top = run_path(str(ROOT / 'ch39-neutral-axial-top.py'))
    boson = run_path(str(ROOT / 'ch39-neutral-bosonic.py'))
    pi = arb.pi()
    g1,g2,v,mz,mp = map(arb,('.467453','.63811','246.22','91.1876','1.22089e19'))
    Q = arb('1.176397e14')
    rho = (v/mp)**2/2
    w,z = rho*g2**2/2,rho*(g2**2+3*g1**2/5)/2
    p0,p1=(mz/mp)**2,(Q/mp)**2
    L_bos = -(2*top['J'](w/p0)+boson['j_neutral'](z/p0))/(4*pi)
    L_bos += boson['scalar_difference'](w/p0,z/p0)/(8*pi)
    a,b = map(arb, ('.1','10'))
    def slope_sign(s):
        return 3*(1+s)*(1/s).log1p()-5
    def f(s):
        return 9*s*(1/s).log1p()-6*s.log1p()
    assert slope_sign(a)>0 and slope_sign(b)<0
    for _ in range(65):
        mid = (a+b)/2
        if slope_sign(mid)>0:
            a=mid
        else:
            assert slope_sign(mid)<0
            b=mid
    s=a.union(b)
    top_bound=f(s)/(144*pi)
    tmax=4*pi*rho
    def ent(m):
        r=m/p1
        return r*(1+(1/(4*r)).log1p())
    high_bos=(2*top['J'](w/p1)+boson['j_neutral'](z/p1))/(4*pi)+(ent(w)+ent(z))/(8*pi)
    high_top=tmax/(4*pi*p1)
    C=(2*w+z)/(24*pi)+(w+z)/(48*pi)+5*tmax/(48*pi)
    finite=25*C/arb('4e-15')
    upper=L_bos+top_bound+high_bos+high_top+finite
    central_t=arb('.8616')**2*rho/p0
    central=L_bos+3*top['J'](central_t)/(8*pi)-top['L'](central_t)/(24*pi)
    for label,x in [('boson low',L_bos),('profile max arg',s),('top upper',top_bound),('boson high upper',high_bos),('top high upper',high_top),('cutoff error',finite),('D upper',upper),('central low D',central)]:
        print(label,x.str(25))
    assert upper < arb('-.051')

    def bosonic_tail(p):
        def entropy(m):
            r=m/p
            return r*(1+(1/(4*r)).log1p())
        return (2*top['J'](w/p)+boson['j_neutral'](z/p))/(4*pi)+(entropy(w)+entropy(z))/(8*pi)
    tail = bosonic_tail((arb('1e13')/mp)**2)
    lower_D_cost=-(L_bos+top_bound+tail)
    needed_inverse=lower_D_cost/C
    min_cutoff=mp*(25/needed_inverse).sqrt()
    for label,x in [('global boson tail',tail),('C bound',C),('inverse moment floor',needed_inverse),('upper bound on lightest cutoff in GeV',min_cutoff)]:
        print(label,x.str(24))
    assert upper < arb('-.0516474') and lower_D_cost-finite > arb('.0516474')
    assert needed_inverse > arb('5.92e32') and min_cutoff < arb('2510')


if __name__ == "__main__":
    verify()

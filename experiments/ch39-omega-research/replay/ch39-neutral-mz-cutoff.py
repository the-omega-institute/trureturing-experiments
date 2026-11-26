"""将最小固有时截止能量取为 M_Z，在原七个较高位置上重求完整端点匹配。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx


def low_weak_column(w, z, t, p0, p1, *, x_over_p0=1):
    # w,z,t 是相对 p0 的质量平方。x_over_p0 指定截止，原算例取一。
    ratio = arb(x_over_p0)
    assert ratio > 0
    def e1(mass):
        return (mass/ratio).expint(1)
    tol = arb("1e-28")
    eta = arb(2)**-60
    def primitive(mass):
        return mass*e1(mass)-ratio*(-mass/ratio).exp()

    ew,et = e1(w),e1(t)
    pw = primitive(w)

    def integrand(s, analytic):
        h = s*(1-s)
        a = s*z
        if analytic and not all(v.real>0 for v in (h,a,a+h,w+h,t+h)):
            return acb("nan")
        charged = w*(ew-e1(w+h))
        neutral = z*(e1(a)-e1(a+h))
        axial = t*(et-e1(t+h))
        vector = 6*h*(e1(h)-e1(t+h))
        scalar = primitive(w+h)-pw-primitive(a+h)+primitive(a)
        return (-(2*charged+neutral)/4+scalar/8+3*axial/8-vector/24)/arb.pi()

    value = acb.integral(integrand,eta,1-eta,abs_tol=tol,rel_tol=tol)
    assert value.is_finite() and value.imag.contains(0)
    # 两类等质量插入 <=h，中性插入 <=1，标量两质量相对零质量核的差各由 h+b 界住，顶矢量核 <=6t。
    # 含系数后的被积函数绝对值 <2，两个尾段的总积分因而至多为 4 eta。
    assert 0<w and 0<z and 0<t<4
    low = value.real+arb(0,4*eta)
    assert low.rad()<arb("1e-14")
    # 高能端用已有的小截止估计，各原始端点的绝对值小于 2x/p。
    return low+arb(0,(2*ratio*p0/p1).upper())


def verify(root):
    original = run_path(str(root/"ch39-full-neutral-free-cutoff.py"))
    gauge = run_path(str(root/"ch39-all-field-signed-crossing.py"))
    gravity = run_path(str(root/"ch39-three-index-free-atoms.py"))
    rotation = run_path(str(root/"ch39-rotation-portal-top-energy.py"))
    gravity_column, _ = gravity["spectra"](root)
    signs,masses = rotation["rotation_signs"](root,2_306_888_100_085)
    data = gauge["joint_system"](root,"1.176397e14",order=12,last=36)
    ctx.prec = 320
    rho = (arb("246.22")/arb("1.22089e19"))**2/2
    portal = arb(".42995401872")
    p0,p1 = data["p0"],data["p1"]
    species, light = original["response_rows"](data,masses,signs,rho,portal)
    pi = arb.pi()
    rule = gauge["angular_rule"](18,60)
    t,w,z = [mass*rho/p0 for _,mass in light]
    weak = low_weak_column(w,z,t,p0,p1)
    assert arb("-.048342336655422") < weak < arb("-.048342336655420")
    print("最低截止对弱三分量的端点差",weak.str(16),flush=True)
    sites = [p0,*map(arb,("4e-15","1e-12","3.5e-11","1e-7","1.23e-6",".5","1"))]
    heat = data["heat"]
    columns,colours = [],[]
    for i,x in enumerate(sites):
        terms = [arb(0) for _ in range(3)]
        for kind,mass,beta in species:
            exposure = gauge["angular_exposure"](kind,mass,x,p0,p1,(p1/p0).log(),rule)/(4*pi)
            for j,b in enumerate(beta):
                terms[j] -= b*exposure
        slope = sum(g*k*x*heat(2,k*rho/x)/4 for g,k in light)
        curvature = -sum(g*k*k*heat(1,k*rho/x)/4 for g,k in light)
        for mass,sign in zip(masses,signs,strict=True):
            shifted = mass+sign*portal*rho
            slope += sign*portal*x*heat(2,shifted/x)
            curvature -= portal**2*heat(1,shifted/x)
        photon,mixed,colour = terms
        columns.append(gravity_column(x)[:3]+[slope/(8*pi**2),curvature/(8*pi**2),
                       photon-8*colour/3,mixed-colour,weak if i==0 else original["weak_column"](x,p0,p1)])
        colours.append(colour)
    a1,a2,a3 = data["initial"]
    target = arb_mat([[n] for n in (1,1,arb("1e-122"),0,arb(".258"),8*a3/3-5*a1/3-a2,a3-a2,0)])
    mat = arb_mat(columns).transpose()
    weights = gravity["solve"](mat,target)
    assert all(not w.contains(0) for w in weights)
    negative = sum(-w for w in weights if w<0)
    colour = a3+sum(w*c for w,c in zip(weights,colours,strict=True))
    assert arb("43.42807") < negative < arb("43.42808")
    assert arb("38.13728") < colour < arb("38.13729")
    determinant = 5*colour**2/3
    assert determinant>0
    print("负权绝对值之和、端点色逆耦合、八个有符号权重",negative.str(12),colour.str(12),
          [w.str(12) for w in weights],flush=True)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)

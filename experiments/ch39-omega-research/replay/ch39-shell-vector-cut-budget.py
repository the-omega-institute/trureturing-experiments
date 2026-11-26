"""最低带荷对与十八个中性超荷重矢量在匹配尺度以下的色散预算。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx

width = run_path(str(Path(__file__).with_name("ch39-shell-yang-mills-vector-width.py")))


def root_upper(x):
    # The physical radicand is nonnegative. Use a pointwise upper bound
    # to include a phase-space endpoint even when an interval straddles it.
    return x.upper().max(arb(0)).sqrt()


def norm(time, space):
    return (time**2+space**2).sqrt()


def spectral_bound(a1, a2, b, s, threshold):
    energy = s.sqrt()
    low = (a2+b)**2
    high = (energy-a1)**2
    span = (s-threshold**2)*(energy-a1+a2+b)/(energy+threshold)
    phase = span**2*root_upper(((energy+a1)**2-low)*(high-(a2-b)**2))/(1024*arb.pi()**2*s*low)

    def kinematics(mass, other_pair):
        maximum = (s+mass**2-other_pair**2)/(2*energy)
        spatial = root_upper(maximum**2-mass**2)
        polarizations = ((arb(0),arb(1),2),(spatial/mass,maximum/mass,1))
        return maximum, spatial, polarizations

    e1,p1,pol1 = kinematics(a1,a2+b)
    e2,p2,pol2 = kinematics(a2,a1+b)
    ez,pz,polz = kinematics(b,a1+a2)
    total = arb(0)
    for first, second, neutral in product(pol1,pol2,polz):
        tz,vz,nz = neutral

        def exchange(spec, rem, ms, mx, es, ps, ex, px):
            ts,vs,_ = spec
            tx,vx,_ = rem
            inner = ((energy-ms)**2-mx**2-b**2)/2
            zz = tz*tx+vz*vx
            left = vs*norm(energy+es,ps)+2*norm(ts,vs)*ps+2*energy*ts
            right = (zz*norm(abs(ez-mx).max(abs(ex-b)),pz+px)
                     +2*norm(tx,vx)*(ex*tz+px*vz).min(root_upper(inner**2/b**2-mx**2))
                     +2*norm(tz,vz)*(ez*tx+pz*vx).min(root_upper(inner**2/mx**2-b**2)))
            longitudinal = (s-ms**2)*abs(b**2-mx**2)*vs*zz/ms**2
            denominator = (b+mx)**2-ms**2
            assert denominator > 0
            return (left*right+longitudinal)/denominator

        t1,v1,n1 = first
        t2,v2,n2 = second
        contact = 2*vz*(t1*t2+v1*v2)+v1*(tz*t2+vz*v2)+v2*(tz*t1+vz*v1)
        amplitude = (contact+exchange(first,second,a1,a2,e1,p1,e2,p2)
                     +exchange(second,first,a2,a1,e2,p2,e1,p1))
        total += n1*n2*nz*amplitude**2
    # For each fixed final polarization the amplitude is a linear functional
    # of the three spatial source components. The preceding bound applies to
    # every unit complex source, hence bounds the squared norm of that
    # functional after squaring. The isotropic source average gives 1/3.
    return total*phase/(6*s)*arb(5)/12


def cut_integral(a,b,steps,second_mass=None):
    other = a if second_mass is None else second_mass
    threshold = a+other+b
    bottom = threshold**2
    total = arb(0)
    for i in range(steps):
        s = bottom+(1-bottom)*arb(i).union(arb(i+1))/steps
        total += (4*spectral_bound(a,other,b,s,threshold)/(s*(1+s))).upper()
    return (1-bottom)*total/steps


def weight_bound(matrix,density,levels,index):
    vector = acb_mat([[v] for v in density])
    factor = arb(1)
    for i,level in enumerate(levels):
        if i != index:
            vector = matrix*vector-vector*level.mid()
            factor *= levels[index]-level.mid()
    return 19*sum(abs(v)**2 for v in vector.entries())/factor**2


def verify():
    ctx.prec = 160
    a, coupling, density, matrix = width["verify"]()
    soft = width["soft"]
    r = soft["geo"]["theta"](804)/soft["geo"]["theta"](12)
    levels = [arb(0)]
    for value in sorted(width["neutral"]["shell"]["shell_spectra"]((0,))[0],key=float):
        if levels[-1].overlaps(value):
            levels[-1] = levels[-1].union(value)
        else:
            levels.append(value)
    full = arb(0)
    for j,level in enumerate(levels[1:],start=1):
        weight = weight_bound(matrix,density,levels,j)
        full += weight*cut_integral(a,level.sqrt()/r,128)
    first = width["neutral"]["invariant_masses"]()[0]
    single = cut_integral(a,first,128)*coupling**2
    print("Lowest invariant mode cut upper per g_D^2",single.upper().str(12))
    print("All eighteen neutral modes cut upper per g_D^2",full.upper().str(12))
    assert single < arb(".00422") and full < arb(".00793")


if __name__ == "__main__":
    verify()

"""三枚相邻破缺谱共用十七节点截止时，从 Z 质量至交汇的全程正定证书。"""
from pathlib import Path
from runpy import run_path
from flint import arb, arb_mat, arb_poly, ctx

SITES = ("1e-50", "1e-19", "1.77827941e-13", "1e-12", "3.16227766e-12", "1e-11",
         "3.5e-11", "1.77827941e-10", "1e-9", "5.623413252e-9", "3.16227766e-8",
         "1.77827941e-7", "2.54e-7", "1e-5", ".0303688", ".5623413251903491", "1")
FIXED = {9: "3830096.77", 6: "-166717.251", 3: "5892.272",
         5: "75053.978", 2: "-571.375"}


def corrected(root):
    def load(short):
        return run_path(str(root/("ch39-"+short+".py")))
    moment, _ = load("three-index-free-atoms")["spectra"](root)
    gauge, neutral = load("all-field-signed-crossing"), load("full-neutral-free-cutoff")
    signs, masses = load("rotation-portal-top-energy")["rotation_signs"](root, 2306888100085)
    data = gauge["joint_system"](root, "1.176397e14", order=4, last=10)
    ctx.prec = 320
    sites = list(map(arb, SITES))
    rho = (arb("246.22")/arb("1.22089e19"))**2/2
    u = arb(".42995401872")
    species, light = neutral["response_rows"](data, masses, signs, rho, u)
    heat, p0, p1 = data["heat"], data["p0"], data["p1"]
    pi, log = arb.pi(), (p1/p0).log()
    rule = gauge["angular_rule"](16, 54)
    cols = []
    for x in sites:
        slope = sum(g*k*x*heat(2, k*rho/x)/4 for g,k in light)
        curvature = -sum(g*k*k*heat(1, k*rho/x)/4 for g,k in light)
        for z,s in zip(masses, signs, strict=True):
            a = (z+s*u*rho)/x
            slope += s*u*x*heat(2,a)
            curvature -= u*u*heat(1,a)
        shifts = [arb(0)]*3
        for kind, z, beta in species:
            val = gauge["angular_exposure"](kind,z,x,p0,p1,log,rule)/(4*pi)
            for j,b in enumerate(beta):
                shifts[j] -= b*val
        qq, qt, colour = shifts
        cols.append(moment(x)+[slope/(8*pi**2), curvature/(8*pi**2), qq-8*colour/3,
                               qt-colour, neutral["weak_column"](x,p0,p1)])
    a1,a2,a3 = data["initial"]
    target = [1,1,arb("1e-122"),0,0,0,0,0,arb(".258"),8*a3/3-(5*a1/3+a2),a3-a2,0]
    # 将弱差中与所有较大截止共有的常数精确地从归一行扣除。
    weak_reference = arb("-.074256577194663365545")
    for col in cols:
        col[-1] -= weak_reference
    target[-1] -= weak_reference
    free = [i for i in range(len(sites)) if i not in FIXED]
    weights = [arb(FIXED[i]) if i in FIXED else arb(0) for i in range(len(sites))]
    b = arb_mat([[v-sum(weights[j]*cols[j][k] for j in FIXED)] for k,v in enumerate(target)])
    mat = arb_mat([cols[j] for j in free]).transpose()
    mid = arb_mat([[v.mid() for v in row] for row in mat.tolist()])
    pre = mid.inv()
    bmid = arb_mat([[v.mid()] for v in b.entries()])
    centers = arb_mat([[v.mid()] for v in mid.solve(bmid).entries()])
    eye = arb_mat([[int(i==j) for i in free] for j in free])
    defect = arb_mat([[abs(v).upper() for v in row] for row in (eye-pre*mat).tolist()])
    residual = [abs(v).upper() for v in (pre*(b-mat*centers)).entries()]
    widths = residual
    for _ in range(12):
        drift = (defect*arb_mat([[x] for x in widths])).entries()
        widths = [(x+y)*arb("1.02")+arb("1e-70") for x,y in zip(residual,drift,strict=True)]
    drift = (defect*arb_mat([[x] for x in widths])).entries()
    assert all(x+y<r for x,y,r in zip(residual,drift,widths,strict=True)), (residual,drift,widths)
    for j,x,r in zip(free, centers.entries(), widths, strict=True):
        weights[j] = x+arb(0,r.upper())
    assert all(not x.contains(0) for x in weights), [x.str(12) for x in weights]
    weak = load("full-neutral-running")["weak_bound"](data,sites,weights,load("neutral-axial-top"),
                                                     load("neutral-bosonic"))
    negative = sum(-w for w in weights if w<0)
    assert arb("6700720.4")<negative<arb("6700720.6") and weak<arb("1.3e-13")
    print("十七权",*[c.str(12) for c in weights], flush=True)
    print("半宽和、总负量、全程弱差界",sum(c.rad() for c in weights).str(8),
          negative.str(13),weak.str(8),flush=True)
    return data, weights, species, weak


def angular_moments(degree):
    # v=4t(1-t)，四个角权均为 A+Bv。返回各自的精确矩。
    raw = [arb(1)]
    for n in range(degree+1):
        raw.append(raw[-1]*(2*n+2)/(2*n+3))
    coeff = [(0,arb(3)/2), (3,-3), (arb(5)/7,arb(3)/7), (arb(9)/11,arb(3)/11)]
    return [[a*raw[n]+b*raw[n+1] for n in range(degree+1)] for a,b in coeff]


def projected(beta):
    return [arb(beta[1]),arb(beta[0])-arb(beta[1]),arb(beta[2])]


def path_certificate(data, weights, species, weak):
    ctx.prec = 1152
    p0 = data["p0"]
    p1 = (arb("1.176397e14")/arb("1.22089e19"))**2
    assert p1.overlaps(data["p1"])
    t0 = p0/p1
    pi = arb.pi()
    sites = list(map(arb,SITES))
    c = [w.mid() for w in weights]
    logspan = (p1/p0).log()
    # 轻场去除共有对数项后只留下正的 Ein 差。按精确归一，所提系数为 1-c[0]。
    # 重场直接展开 E1(z+hp) 的差。所有近似都分别控制上下两端。
    heavy_order, light_order = 56, 510
    angular = angular_moments(light_order)
    heavy = [[arb_poly([]) for _ in range(4)] for j in range(3)]
    light = [[arb(0)]*(light_order+1) for j in range(3)]
    err = [arb(0)]*3
    medium = [[arb(0)]*3 for _ in range(4)]
    m = p1/8
    ein_size = [(p1/(4*x)).min(1+(1+p1/(4*x)).log()) for x in sites]
    for kind,z,beta in species:
        beta = projected(beta)
        if z < arb("1e-25"):
            # 重组前保留精确权的归一，使用轻质量单节点 2z/x 误差。
            unc = abs(1-weights[0])*logspan+abs(weights[0])*12*sites[0]/p0
            unc += sum(abs(w)*2*z/x+w.rad()*size
                       for w,x,size in zip(weights[1:],sites[1:],ein_size[1:],strict=True))
            for j,b in enumerate(beta):
                err[j] += (abs(b)*unc+b.rad()*sum(abs(k)*size for size,k in
                           zip(ein_size[1:],c[1:],strict=True)))/(4*pi)
                medium[kind][j] += b.mid()*c[1]/(4*pi)
            for x,k in zip(sites[2:],c[2:],strict=True):
                a = p1/(4*x)
                term = arb(1)
                for n in range(1,light_order+1):
                    term *= -a/n
                    val = -k*term*angular[kind][n]/(n*4*pi)
                    for j,b in enumerate(beta):
                        light[j][n] += b.mid()*val
                tail = abs(term)*a/(light_order+1)**2
                assert a<light_order+1
                for j,b in enumerate(beta):
                    err[j] += abs(k*b)*tail*2/(4*pi)
            continue
        assert z>m
        # s=2vq-1，q=p/p1。在 |s|=2 上 z+m+ms 的实部至少 z-m。
        # 若最小节点的指数极小，直接在物理实区间抑制整项。
        for j,b in enumerate(beta):
            err[j] += abs(b)*sum(w.rad()*(-z/x).exp() for w,x in
                                 zip(weights,sites,strict=True))*p1/(16*pi*z)
        hc = [arb(0)]*(heavy_order+1)
        tail = arb(0)
        for x,k in zip(sites,c,strict=True):
            if z/x>150:
                tail += abs(k)*logspan*(-z/x).exp()
                continue
            tail += abs(k)*4*m/(z-m)*arb(2)**(-heavy_order)
            a,d = m/x, m/(z+m)
            source = (-((z+m)/x)).exp()
            val = source
            for n in range(1,heavy_order+1):
                hc[n] -= k*d*val/n
                source *= -a/n
                val = source-d*val
        pol = arb_poly(hc)(arb_poly([-1,2]))
        for j,b in enumerate(beta):
            heavy[j][kind] += pol*(b/(4*pi))
            err[j] += abs(b)*tail/(4*pi)
    poly = [arb_poly(q) for q in light]
    for j in range(3):
        for k in range(4):
            poly[j] += arb_poly([v*angular[k][n] for n,v in enumerate(heavy[j][k])])
    a1,a2,a3 = data["initial"]
    initial = (a2,5*a1/3,a3)
    baseline = [v-f(t0) for v,f in zip(initial,poly,strict=True)]
    # 第二个节点跨越多个能标；对 Ein 角平均采用正权单调求积。
    bins = 512
    mesh = [arb(i)/(2*bins) for i in range(bins+1)]
    h = [t*(1-t) for t in mesh]
    def angle_cdf(t):
        v = 4*(t*t/2-t**3/3)
        return (3*v,6*t-6*v,(10*t+6*v)/7,(18*t+6*v)/11)
    cdf = [angle_cdf(t) for t in mesh]
    cellweights = [[v-u for u,v in zip(a,b,strict=True)] for a,b in zip(cdf[:-1],cdf[1:],strict=True)]
    def ein(y):
        if y.is_zero():
            return arb(0)
        if y>100:
            return y.log()+arb.const_euler()+arb(0,arb("4e-46"))
        return y.log()+arb.const_euler()+y.expint(1)
    memo = {}
    def medium_values(q):
        values = [ein(q*p1*s/sites[1]) for s in h]
        return [sum(cell[k]*l for cell,l in zip(cellweights,values[:-1],strict=True)).union(
                sum(cell[k]*r for cell,r in zip(cellweights,values[1:],strict=True))) for k in range(4)]
    zero = medium_values(t0)
    print("全程多项式近似误差",*[v.str(8) for v in err],flush=True)
    def value(num, depth):
        key = arb(num)*arb(2)**(-depth)
        # 用整数配对复用同一个二进端点。
        reduced = (num,depth)
        if num==0:
            reduced = (0,0)
        else:
            while reduced[1] and reduced[0]%2==0:
                reduced = (reduced[0]//2,reduced[1]-1)
        if reduced not in memo:
            memo[reduced] = medium_values(key)
        return memo[reduced]
    floors = [arb(1000)]*3
    stack = [(0,0)]
    leaves = 0
    while stack:
        num,dep = stack.pop()
        left,right = value(num,dep),value(num+1,dep)
        half = arb(2)**(-dep-1)
        mid = (2*num+1)*half
        bounds=[]
        for j,f in enumerate(poly):
            localized = f(arb_poly([mid,half]))
            val = localized[0]-sum(abs(t).upper() for t in localized.coeffs()[1:])
            val += baseline[j]-err[j]+sum(medium[k][j]*(left[k].union(right[k])-zero[k]) for k in range(4))
            bounds.append(val.lower())
        if min(bounds)>1:
            floors = [min(a,b) for a,b in zip(floors,bounds,strict=True)]
            leaves+=1
        else:
            assert dep<80,(num, dep,[a.str(9) for a in bounds], [f(mid).str(9) for f in poly])
            stack.extend(((num*2+1,dep+1),(num*2,dep+1)))
    b,y,colour = floors
    assert b>arb("17.40") and y>arb("3.47") and colour>arb("1.389")
    assert b*y/(b+y)>weak
    print("覆盖闭段数及 QT3、QQ−QT3、色的严格下界",leaves,*[v.str(11) for v in floors],flush=True)


if __name__=="__main__":
    path_certificate(*corrected(Path(__file__).resolve().parent))

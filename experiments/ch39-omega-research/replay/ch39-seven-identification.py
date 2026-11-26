"""七条件交点的区间雅可比、盒内唯一性及轻耦合灵敏度。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, arb_mat, ctx


def dot(a, b):
    return sum(x*y for x, y in zip(a, b, strict=True))


def check_bounds(value, lower, upper):
    assert arb(lower) < value < arb(upper), (value, lower, upper)


def projected(beta):
    b1, b2, b3 = map(arb, beta)
    return (5*b1/3+b2-8*b3/3, b2-b3)


def energy_derivative(kind, z, x, p, rule):
    # dR/d log M 的积分。实轴上被积核在 [0,2 exp(-z/x)] 内。
    upper = 2*(-z/x).exp()
    if upper < arb("1e-16"):
        return arb(0).union(upper)
    nodes, error = rule
    total = arb(0)
    for h, weights in nodes:
        ph = p*h
        a = z+ph
        ratio = arb(1) if z.is_zero() else ph/a
        total += weights[kind]*ratio*(-a/x).exp()
    return 2*(total+arb(0, error.upper()))


def weak_mass_derivatives(zw, p0, p1):
    # 无限截止下 z*(-dR/dz) 的公共项。各正角权的单端点矩随 z/p 增大。
    def endpoint(s):
        a = 1+4*s
        plain = 4*s/a.sqrt()*(1/a.sqrt()).atanh()
        angular = a*plain-4*s
        return 3*(4*plain-angular)/11, 3*angular

    low = endpoint(zw/p0)
    ratio = zw/p1
    lower = endpoint(ratio.lower())
    upper = endpoint(ratio.upper())
    return [start-lo.union(hi) for start, lo, hi in zip(low,lower,upper,strict=True)]


def verify(root):
    ctx.prec = 192
    pi = arb.pi()
    energy = arb("1.1763969e14").union(arb("1.17639745e14"))
    portal = arb(".42995401871").union(arb(".42995401873"))
    g2, yt, g1 = map(arb, (".63811", ".8616", ".467453"))
    six = run_path(str(root/"ch39-six-broken-conditions.py"))
    system = six["broken_system"](root, energy, portal)
    data = system["original"]
    weights = system["weights_at"](portal)[0]
    rho, zw = system["rho"], system["weak_mass"]
    p0, p1 = data["p0"], data["p1"]
    heat = data["heat"]
    masses, signs = system["masses"], system["signs"]
    zero = arb(0)
    # z=z0+k*rho；最后三列依次为 k 对 u、g2、yt 的偏导。
    moving = [(-12, arb(1), zero, yt**2, (0,0,2*yt)),
              (6, arb(-1), zero, g2**2/2, (0,g2,0)),
              (3,-arb(1)/2, zero, (g2**2+3*g1**2/5)/2, (0,g2,0))]
    moving.extend((4,arb(2)/3,z,s*portal,(s,0,0)) for z,s in zip(masses,signs,strict=True))

    fermion, scalar, _, gauge = range(4)
    w_gauge = (-arb(22)/3,-arb(22)/3)
    w_scalar = (arb(1)/3,arb(1)/6)
    species = [(fermion,zero,projected(data["species"][0][2])),
               (gauge,zero,(arb(88)/3,arb(11))),
               (gauge,zw,w_gauge),(scalar,zw,w_scalar)]
    j = 0
    for kind, z, beta in data["species"][3:]:
        if kind == scalar:
            assert z.overlaps(masses[j])
            z += signs[j]*portal*rho
            j += 1
        species.append((kind,z,projected(beta)))
    assert j == 16
    angular = run_path(str(root/"ch39-all-field-signed-crossing.py"))["angular_rule"](14,39)
    weak = weak_mass_derivatives(zw,p0,p1)
    heavy_bound = sum(rho/z.union(z+s*portal*rho).lower() for z,s in zip(masses,signs,strict=True))
    log_moment = run_path(str(root/"ch39-mixed-neutral.py"))["mixed_log_moment"]
    r1 = zw/p1
    light_shift = (log_moment(r1.lower()).union(log_moment(r1.upper()))-log_moment(zw/p0))/(4*pi)

    matrix, slope_cols = system["columns_at"](portal)
    assert not matrix.det().contains(0)
    c_cols, partials = [], [[],[],[],[]]
    for x, sym in zip(data["sites"],data["cols"],strict=True):
        # 每列依次是归一、K、W、B、U''、C、U' 的单截止贡献。
        deriv = [[zero]*7 for _ in range(4)]
        for g, alpha, z0, kappa, change in moving:
            z = z0+kappa*rho
            assert z > 0
            e1, e2 = heat(1,z/x), heat(2,z/x)
            for j, dk in enumerate(change,1):
                if dk == 0:
                    continue
                deriv[j][1] -= alpha*rho*dk*e1/(2*pi)
                deriv[j][2] += g*x*rho*dk*e2/(4*pi)
                deriv[j][6] += g*dk*(x*e2-kappa*rho*e1)/(32*pi*pi)
                deriv[j][4] += g*dk*(-2*kappa*e1+kappa**2*rho*(-z/x).exp()/z)/(32*pi*pi)
        for kind, z, coeff in species:
            de = energy_derivative(kind,z,x,p1,angular)/(4*pi)
            deriv[0][3] -= coeff[0]*de
            deriv[0][5] -= coeff[1]*de
        for row, gauge_beta, scalar_beta in zip((3,5),w_gauge,w_scalar,strict=True):
            # 标量门户的 B、C 荷恰与相应带电 Goldstone 的荷相等。
            deriv[1][row] = arb(0,(scalar_beta*heavy_bound/(4*pi)).upper())
            weak_terms = [v+(-zw/x).union(0) for v in weak]
            deriv[2][row] = 2*(gauge_beta*weak_terms[0]+scalar_beta*weak_terms[1])/(g2*4*pi)
        mixed_error = ((arb(43)/6)*zw/x+abs(portal)*heavy_bound/6)/(4*pi)
        mixed = sym[4]-sym[5]+light_shift+arb(0,mixed_error.upper())
        c_cols.append(mixed)
        for ds, dest in zip(deriv,partials,strict=True):
            dest.append(ds)

    # 归一化给 sum dc=0；减去一个共同常数可收紧 C 的导数。
    pivot = c_cols[2].mid()
    c_cols = [v-pivot for v in c_cols]
    derivatives, weight_derivatives = [], []
    for part in partials:
        rhs = arb_mat([[-dot(weights,[col[i] for col in part])] for i in range(5)])
        dc = matrix.solve(rhs).entries()
        weight_derivatives.append(dc)
        derivatives.append([dot(dc,c_cols)+dot(weights,[col[5] for col in part]),
                            dot(dc,slope_cols)+dot(weights,[col[6] for col in part])])
    jac = arb_mat(derivatives[:2]).transpose()
    params = arb_mat(derivatives[2:]).transpose()
    print("闭盒雅可比和行列式",jac, jac.det().str(14),flush=True)
    for value, bounds in zip(jac.entries(),
            (("-3.466","-3.464"),("23.81","23.83"),("-9.873e-11","-9.871e-11"),
             ("-.018542093","-.018542091")),strict=True):
        check_bounds(value,*bounds)
    check_bounds(jac.det(),".06424",".06426")
    sensitivity = jac.solve(-params)
    print("解对 g2 与 yt 的绝对输入导数",sensitivity,flush=True)
    bounds = [(("-1.355","-1.353"),("-.229041158","-.229041156")),
              ((".67",".73"),("1.113531853","1.113531855"))]
    for i, coupling in enumerate((g2,yt)):
        for j in range(2):
            check_bounds(coupling*sensitivity[j,i],*bounds[i][j])
        print("对数输入导数（log M、u、五权）",*(z.str(12) for z in
              [coupling*sensitivity[j,i] for j in range(2)]+
              [coupling*(weight_derivatives[2+i][k]+weight_derivatives[0][k]*sensitivity[0,i]
                         +weight_derivatives[1][k]*sensitivity[1,i]) for k in range(5)]),flush=True)

if __name__=="__main__":
    root = Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
    del sys.argv[1:]
    verify(root)

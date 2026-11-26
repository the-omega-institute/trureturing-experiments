"""弱中心互逆刚度在任意固定标量边相位上的统一一圈交点。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx, fmpz_mat

center = run_path(str(Path(__file__).with_name("ch39-native-center-parity-crossing.py")))


def dot(v, w):
    return sum(a*b for a, b in zip(v, w))


def verify():
    ctx.prec = 192
    shell, mz, pi = center["shell"], center["mz"], arb.pi()
    soft = shell["soft"]
    height = acb.zeta_zero(11431).imag
    base, neutral, _ = soft["vector_pencil"](height)
    graph = base + neutral
    pole = soft["pole"]
    keep = [i for i in range(19) if i != pole]
    assert all(graph[pole, i] == -1 for i in keep)
    assert sum(graph[i, i] for i in range(19)) == 270
    pinned = fmpz_mat([[int(graph[i, j].real.unique_fmpz()) for j in keep] for i in keep])
    w = 1 - arb(23)/25*arb(2)**-20/(height**2 + soft["point"]["study"]["offset"])
    poles = [w*z.real for z, k in pinned.charpoly().complex_roots() for _ in range(k)]
    weights = shell["shell_spectra"]((0,))[0]
    radius = soft["geo"]["theta"](804)/soft["geo"]["theta"](12)
    r2 = radius**2
    Lmin, Lmax = [(arb(mu)/mz).log() for mu in ("1e13", "1e17")]
    reference = r2*(-2*Lmin).exp()
    umin, umax = arb("2.42"), arb("3.21")
    # The retained 18 eigenvalues are >= w by the spanning star. Their
    # mean is <=15 by the trace, and all eigenvalues are <=36 by degree.
    assert all(v > 0 for v in poles)
    fer = lambda mass: center["sheet"](mass, r2, False)
    sca = lambda mass: center["scalar"](mass, r2)
    vec = lambda mass: center["sheet"](mass, r2, True)

    def exposure(lo, hi, kernel):
        # Jensen bounds the sum from below. Principal-submatrix interlacing,
        # heat-kernel diamagnetism and the edge minimum w bound it from above.
        ref = 18*reference/(8*w*lo)
        return (18*kernel(15*hi) - ref).union(sum(kernel(lo*v) for v in poles))

    def box(lo, hi):
        odd = exposure(1/hi, 1/lo, fer)
        even = exposure(lo, hi, fer)
        scalar = exposure(1/hi, 1/lo, sca)
        vector = (sum(vec(hi*v) for v in weights) - 18*reference/(8*w*lo)).union(
            sum(vec(lo*v) for v in weights))
        return [odd, even, even, odd, even, scalar, vector]

    em, weak = arb("127.930"), arb(".23129")
    initial = [3*(1-weak)*em/5, weak*em, 1/arb(".1180")]
    rows = [[-3*v/(2*pi) for v in row] for row in center["weyl"]]
    rows += [[-v/(2*pi) for v in center["higgs"]], [arb(0), 7/(2*pi), 21/(4*pi)]]
    P = [arb(v)/218 for v in (115, -333, 218)]
    Y = [1, -1, 0]
    pp, yy = [dot(P, row) for row in rows], [dot(Y, row) for row in rows]
    # Odd exposures increase and even ones decrease with log u.
    signs = (1, -1, -1, 1, -1, 1, -1)
    assert all(sign*c > 0 for c, sign in zip(pp, signs))
    assert dot(P, initial) + dot(pp, box(umin, umin)) < 0
    assert dot(P, initial) + dot(pp, box(umax, umax)) > 0

    # Y at P=0 can equally be computed using Y-kP. The ratios remove one
    # interval term apiece. Intersect all resulting enclosures.
    eliminators = [y/p for y, p in zip(yy, pp)]
    beta = center["beta"]
    slope = (beta[0]-beta[1])/(2*pi)
    log_lo, log_hi = [], []
    common_lo, common_hi = [], []
    for i in range(96):
        lo, hi = [umin + (umax-umin)*j/96 for j in (i, i + 1)]
        exposures = box(lo, hi)
        if not (dot(P, initial) + dot(pp, exposures)).contains(0):
            continue
        candidates = []
        for k in eliminators:
            projection = [y-k*p for y, p in zip(Y, P)]
            coeff = [dot(projection, row) for row in rows]
            candidates.append((dot(projection, initial) + dot(coeff, exposures))/slope)
        low = max(v.lower() for v in candidates)
        high = min(v.upper() for v in candidates)
        if high < low:
            continue
        log_lo.append(low)
        log_hi.append(high)
        common = initial[0] + dot([row[0] for row in rows], exposures) - beta[0]*low.union(high)/(2*pi)
        common_lo.append(common.lower())
        common_hi.append(common.upper())
    fitted = min(log_lo).union(max(log_hi))
    mu = mz*fitted.exp()
    common = min(common_lo).union(max(common_hi))
    assert arb("4.0826e13") < mu < arb("7.8786e13")
    assert arb("37.4757") < common < arb("38.5175")

    # At each fixed L the sign of every summand gives a unique P root.
    # Bound its variation in L before solving the second crossing equation.
    ref_fraction = reference/(4*w*min(umin, 1/umax))
    score = (3*r2/16)/(36/umin + 3*r2/16)
    p_x = pp[0]*arb(11)/32*(score - ref_fraction)
    p_L = 18*sum(abs(c) for c in pp)*ref_fraction
    y_L = 18*sum(abs(c) for c in yy)*ref_fraction
    y_x = 9*sum(abs(c) for c in yy)
    assert p_x > 0 and y_L + y_x*p_L/p_x < slope/100

    def hypercharge_at(factor):
        # The fitted masses are now fixed. The two momenta are mu and factor*mu,
        # so the M_Z subtraction cancels exactly from this difference.
        def step(scale_lo, scale_hi, kernel):
            def difference(mass):
                return kernel(mass, r2*factor**2) - kernel(mass, r2)
            return (18*difference(15*scale_hi)).union(
                sum(difference(scale_lo*v) for v in poles))
        fermion = lambda mass, momentum: center["sheet"](mass, momentum, False)
        odd = step(1/umax, 1/umin, fermion)
        even = step(umin, umax, fermion)
        scalar = step(1/umax, 1/umin, center["scalar"])
        assert (3*(center["weyl"][0][0] + center["weyl"][3][0]) - arb(4)/5).contains(0)
        assert (3*sum(center["weyl"][j][0] for j in (1, 2, 4)) - arb(16)/5).contains(0)
        loss = (4*odd/5 + 16*even/5 + scalar/10 + beta[0]*arb(factor).log())/(2*pi)
        return common-loss
    assert hypercharge_at(30) > 0 > hypercharge_at(45)
    print("Uniform stiffness interval", umin.str(8), umax.str(8))
    print("Uniform matching scale bounds / GeV", mu.lower().str(12), mu.upper().str(12))
    print("Uniform common inverse bounds", common.lower().str(12), common.upper().str(12))
    print("Uniform bracket for the first hypercharge inverse zero / matching scale", 30, 45)


if __name__ == "__main__":
    verify()

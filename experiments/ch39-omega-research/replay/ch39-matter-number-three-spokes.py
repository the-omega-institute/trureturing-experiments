"""零裸质量的 75 维物质在三条辐边上携带共同粒子数相位。"""
from collections import Counter
from fractions import Fraction
from functools import cache
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, acb_poly, arb, arb_poly, ctx, fmpq, fmpq_poly, fmpz_mat

crossing = run_path(str(Path(__file__).with_name("ch39-unit-seventy-five-crossing.py")))
smooth, soft = crossing["smooth"], crossing["soft"]



def spoke_polynomials(matrix, pole, orbit, weight, relative):
    def polynomial(phase):
        m = acb_mat(matrix)
        for j, factor in zip(orbit, (relative, relative.conjugate(), 1)):
            m[pole,pole] += weight-1
            m[j,j] += weight-1
            m[pole,j] = matrix[pole,j]*weight*phase*factor
            m[j,pole] = m[pole,j].conjugate()
        coefficients = (-m).charpoly().coeffs()
        assert all(z.imag.contains(0) for z in coefficients)
        return arb_poly([z.real for z in coefficients])
    plus, minus, side = (polynomial(z) for z in (1, -1, acb(0,1)))
    constant = (plus+minus)*arb(".5")
    return constant, (plus-minus)*arb(".5"), side-constant


def verify_vacuum(matrices, pole, orbit, weight):
    polynomials = {}
    for q, matrix in matrices.items():
        constant, cosine, sine = spoke_polynomials(matrix, pole, orbit, weight, acb(1))
        # A phase-dependent permutation uses an oriented cycle on >=3
        # vertices; s^17, s^18 and s^19 are exactly phase independent.
        for poly in (cosine, sine):
            assert all(poly[k].contains(0) for k in (17,18,19))
        cosine = arb_poly(cosine.coeffs()[:17])
        sine = arb_poly(sine.coeffs()[:17])
        if q == 0:
            assert len(cosine.coeffs()) == 17 and all(t < 0 for t in cosine.coeffs())
        else:
            margin = -constant*cosine-(cosine**2+sine**2)
            assert len(margin.coeffs()) == 36 and all(t > 0 for t in margin.coeffs())
            curvature = margin**2-sine**2*(constant**2-cosine**2-sine**2)
            assert len(curvature.coeffs()) == 71 and all(t > 0 for t in curvature.coeffs())
        polynomials[q] = (constant, cosine, sine)
    print("Strict free determinant monotonicity and strict convexity of the negative logarithm in cos(angle)")
    return polynomials


def integrate_force(numerator, denominator, tolerance="1e-11"):
    # Integrate s*N(s)/D(s) in s=16*t and s=16/t. The tail singularity
    # cancels since deg(D)-deg(N)>=3.
    n, d = numerator.degree(), denominator.degree()
    assert n <= d-3
    def mapped(poly, reverse):
        coeff = [v*16**i for i,v in enumerate(poly.coeffs())]
        return acb_poly(coeff[::-1] if reverse else coeff)
    result = arb(0)
    for reverse in (False, True):
        top, bottom = mapped(numerator, reverse), mapped(denominator, reverse)
        def integrand(t, analytic):
            return 256*(t**(d-n-3) if reverse else t)*top(t)/bottom(t)
        val = acb.integral(integrand, 0, 1, rel_tol=arb(tolerance), abs_tol=arb(tolerance)/10)
        assert val.imag.contains(0) and val.is_finite()
        result += val.real
    return result


def wilson_force(polynomials, cosine):
    total = arb(0)
    for q, multiplicity in ((0,33), (5,18), (10,3)):
        a,b,d = polynomials[q]
        if q == 0:
            numerator, denominator = -b, a+b*cosine
        else:
            norm = b**2+d**2
            numerator = -2*(a*b+norm*cosine)
            denominator = a**2-d**2+2*a*b*cosine+norm*cosine**2
        total += multiplicity*integrate_force(numerator,denominator)
    return total/(8*arb.pi()**2)


def relative_curvature(matrices, pole, orbit, weight, c):
    # The cyclic permutation fixes the common configuration; its derivative
    # in the (1,-1,0) direction is zero separately in every charge tower.
    soft["star"]["charged"]["geometry"]()
    sine = (1-c*c).sqrt()
    total = arb(0)
    for q, matrix in matrices.items():
        values = []
        for relative in (1, -1, acb(0,1), acb(0,-1)):
            a,b,d = spoke_polynomials(matrix,pole,orbit,weight,acb(relative))
            values.append(a+b*c+d*sine)
        start, opposite, left, right = values
        # Angular Fourier degree <=2; these four values determine -P''(0).
        numerator = (start-opposite)*arb(".5")+(start+opposite-left-right)
        assert all(numerator[k].contains(0) for k in (17,18,19))
        numerator = arb_poly(numerator.coeffs()[:17])
        total += {0:33,5:18,10:3}[abs(q)]*integrate_force(numerator,start,"2e-5")
    return total/(8*arb.pi()**2)


def other_spoke_forces(matrices, pole, orbit, weight, c, neighbors):
    sine = (1-c*c).sqrt()
    def determinant(matrix):
        a,b,d = spoke_polynomials(matrix,pole,orbit,weight,acb(1))
        return a+b*c+d*sine
    denominators = {q:determinant(matrix) for q,matrix in matrices.items()}
    rows = []
    for extra in soft["star"]["charged"]["orbits"][1:]:
        free = arb(0)
        for q,matrix in matrices.items():
            moved = []
            for sign in (1,-1):
                m = acb_mat(matrix)
                for j in extra:
                    assert m[pole,j] == -1
                    m[pole,j] = acb(0,-sign)
                    m[j,pole] = m[pole,j].conjugate()
                moved.append(determinant(m))
            numerator = (moved[1]-moved[0])*arb(".5")
            assert all(numerator[k].contains(0) for k in (17,18,19))
            numerator = arb_poly(numerator.coeffs()[:17])
            free += {0:33,5:18,10:3}[abs(q)]*integrate_force(numerator,denominators[q],"2e-5")
        triangles = sum(k in neighbors[j] for k in orbit for j in extra)
        squares = sum(len(neighbors[j]&neighbors[k]) for k in orbit for j in extra)
        rows.append((triangles,squares,free/(8*arb.pi()**2*sine)))
    return rows


def verify():
    ctx.prec = 224
    h = acb.zeta_zero(11431).imag
    base, neutral, _ = soft["vector_pencil"](h)
    pole, orbit = soft["pole"], soft["star"]["orbit"]
    assert len(orbit) == 3
    matter = crossing["matter"]
    signed = {q: base+matter["nonspokes"](q,h) for q in (-10,-5,5,10)}
    ordinary = fmpz_mat([[int(v.real.unique_fmpz()) for v in row]
                         for row in (base+neutral).tolist()])
    p_plus = fmpq_poly(ordinary.charpoly())
    for j in orbit:
        assert ordinary[pole,j] == ordinary[j,pole] == -1
        ordinary[pole,j] = ordinary[j,pole] = 1
    p_minus = fmpq_poly(ordinary.charpoly())
    change = arb(23)/25*arb(2)**-20/(h*h+soft["point"]["study"]["offset"])
    polynomials = verify_vacuum({0:base+neutral, 5:signed[5], 10:signed[10]}, pole, orbit, 1-change)
    weyl = arb(0).union(4*change)
    geo = soft["geo"]
    r2 = (geo["theta"](804)/geo["theta"](12))**2
    vector_mass = crossing["shell_spectra"]((0,5))
    v0, v5 = (sum(smooth["sheet"](v,r2,True) for v in vector_mass[q]) for q in (0,5))
    mz, em, weak = arb("91.1876"), arb("127.930"), arb(".23129")
    initial = (3*(1-weak)*em/5, weak*em, 1/arb(".1180"))
    pi, beta = arb.pi(), smooth["beta"]
    background = [v+7*(n*v0+k*v5)/(4*pi) for v,n,k in zip(initial,(0,2,3),(5,3,2))]
    weights = crossing["index_table"]()

    @cache
    def eigen(cosine):
        c = fmpq(cosine.numerator, cosine.denominator)
        phase = acb(arb(c), (1-arb(c)**2).sqrt())
        # Every determinant monomial meets the pole row and column at most
        # once each. Reality and conjugation leave an affine function of cos.
        polynomial = ((1+c)*p_plus+(1-c)*p_minus)/2
        roots = polynomial.complex_roots()
        assert all(z.imag.contains(0) for z,k in roots)
        result = {0:[z.real-weyl for z,k in roots for _ in range(k)]}
        for q, matrix in signed.items():
            m = acb_mat(matrix)
            for j in orbit:
                m[pole,j], m[j,pole] = -phase, -phase.conjugate()
            result[q] = [v-weyl for v in smooth["masses"](m)]
        assert all(len(row) == 19 for row in result.values())
        return result

    def intercept(cosine, allowance=0):
        def read(row):
            # Weyl errors can be inserted into the mass before taking sheet.
            assert all(v-allowance > arb(".0001") for v in row)
            return sum(smooth["sheet"](v+arb(0,allowance), r2, False) for v in row)
        d = {q: read(row) for q,row in eigen(cosine).items()}
        return [v-sum(row[i]*(d[q]+d[-q]) for q,row in weights.items())/(3*pi)
                for i,v in enumerate(background)]

    def projection(a):
        return (beta[1]-beta[2])*(a[0]-a[2])-(beta[0]-beta[2])*(a[1]-a[2])

    a, b = Fraction(9,10), Fraction(23,25)
    sign = projection(intercept(a))
    assert sign*projection(intercept(b)) < 0
    for _ in range(30):
        mid = (a+b)/2
        side = projection(intercept(mid))
        assert not side.contains(0)
        if side*sign > 0:
            a = mid
        else:
            b = mid
    # Use rational endpoints and a centered operator norm enclosure for all
    # intermediate cosines. The twisted off-diagonal star has norm sqrt(3).
    a -= Fraction(1,10**9)
    b += Fraction(1,10**9)
    center = (a+b)/2
    def ball(v):
        return arb(fmpq(v.numerator,v.denominator))
    motion = arb(3).sqrt()*ball(b-a)/(2*(1-ball(b)**2).sqrt())
    whole = intercept(center, motion.upper())
    L = 2*pi*(whole[0]-whole[2])/(beta[0]-beta[2])
    logs = L+arb(0,arb("1e-5"))
    ref2 = r2*(-2*logs).exp()
    err = 19*ref2/(8*arb(".0001"))*(35/(4*pi)+50/(3*pi))
    error = arb(0,err.upper())
    assert projection([v+error for v in intercept(a)])*projection([v+error for v in intercept(b)]) < 0
    for t,sgn in ((logs.lower(),1),(logs.upper(),-1)):
        assert sgn*(whole[0]-whole[2]-(beta[0]-beta[2])*t/(2*pi)+2*error) > 0
    inverse = whole[2]-beta[2]*logs/(2*pi)+error
    angle = ball(a).union(ball(b)).acos()
    value = mz*logs.exp()
    assert inverse > 17 and arb("1e16") < value < arb("2e16")
    sites = geo["sites"]
    distances = {sum((u-v)**2 for u,v in zip(x,y)) for i,x in enumerate(sites) for y in sites[:i]}
    assert len(distances) == 30
    radius_sq = sum(u*u for u in sites[0])
    assert all(sum(u*u for u in x) == radius_sq for x in sites)
    window_ref2 = r2*(mz/arb("1e13"))**2
    window_err = arb(0,(19*window_ref2/(8*arb(".0001"))*(35/(4*pi)+50/(3*pi))).upper())
    gaps = []
    for d in sorted(distances):
        direct = intercept(1-Fraction(d,2*radius_sq))
        p = projection([v+window_err for v in direct])
        # The positive coefficients of this zero-sum projection total bY-bC.
        gaps.append(abs(p).lower()/(beta[0]-beta[2]))
    assert min(gaps) > arb(".05222")
    at_pi = projection([v+window_err for v in intercept(Fraction(-1))])
    vacuum_gap = abs(at_pi)/(beta[0]-beta[2])
    assert vacuum_gap > arb("2.51")
    print("Spread floor at the free vacuum", vacuum_gap.lower().str(11))
    assert arb(".43546463") < angle < arb(".43546465")
    assert arb("1.3676e16") < value < arb("1.3678e16")
    assert arb("17.1595") < inverse < arb("17.1597")
    print("Three-spoke matter phase (rad), scale GeV, common inverse", angle.str(12),value.str(12),inverse.str(10))
    print("Global spread floor for all native single pair angles", min(gaps).lower().str(11))

    force = wilson_force(polynomials,ball(a)).union(wilson_force(polynomials,ball(b)))
    edges = {(i,j) for i,j,_,_ in geo["edges"]}
    assert all(tuple(sorted((pole,j))) in edges for j in range(len(sites)) if j != pole)
    faces = sum(pole not in edge and len(set(edge).intersection(orbit)) == 1 for edge in edges)
    print("Flux carrying triangles, total Wilson coefficient, coefficient per face", faces,
          force.str(12), (force/faces).str(12))
    assert faces == 27 and arb(".10231448") < force/faces < arb(".10231450")
    assert all(not set(edge).issubset(orbit) for edge in edges)
    cbox = ball(a).union(ball(b))
    relative = relative_curvature({0:base+neutral, **signed},pole,orbit,1-change,cbox)
    total_relative = relative + 2*cbox*force/3
    print("Curvature for angular changes (1,-1,0), free and with triangles", relative.str(11),total_relative.str(11))
    assert -arb(".1076") < total_relative < -arb(".1074")
    # One unoriented simple square through the pole is fixed by its two
    # neighbours there and their common intermediate site away from the pole.
    others = set(range(len(sites)))-{pole}
    neighbors = {j:{k for k in others if tuple(sorted((j,k))) in edges} for j in others}
    counts = Counter()
    for j,k in combinations(others,2):
        counts[frozenset((j,k)).intersection(orbit)] += len(neighbors[j] & neighbors[k])
    assert all(counts[frozenset((j,))] == 114 for j in orbit)
    assert all(counts[frozenset(jk)] == 5 for jk in combinations(orbit,2))
    single = faces+3*114
    both_weight = force/single
    both_relative = total_relative+30*both_weight
    assert arb(".1170") < both_relative < arb(".1173")
    print("Triangles and squares: coefficient per loop, transverse curvature",
          both_weight.str(11),both_relative.str(11))
    rows = other_spoke_forces({0:base+neutral, **signed},pole,orbit,1-change,cbox,neighbors)
    for n3,n4,derivative in rows:
        print("Extra radial orbit: triangle and square count, free angular force / sine",
              n3,n4,derivative.str(10))
    orbits = soft["star"]["charged"]["orbits"]
    swap = {9:10, 10:9, 13:16, 16:13, 15:17, 17:15}
    idx = lambda j: swap.get(j,j)
    assert {tuple(sorted((idx(j),idx(k)))) for j,k in edges} == edges
    assert idx(pole) == pole and all(idx(j) == j for j in orbit)
    assert {idx(j) for j in orbits[3]} == set(orbits[4])
    assert rows[2][2]-rows[3][2] > arb(".0008")
    # The differential Ward identity is exact: a common shift of every spoke
    # is a vertex gauge transformation at the pole.
    assert (sum(row[2] for row in rows)-force).contains(0)

    lengths = {(i,j): d for i,j,d,_ in geo["edges"]}
    def length(i,j):
        return lengths[tuple(sorted((i,j)))]
    shapes = []
    for extra in orbits[1:]:
        shapes.append(Counter(tuple(sorted((length(pole,i),length(pole,j),length(i,j))))
                              for i in orbit for j in extra if j in neighbors[i]))
    unique = [(32,544,544),(68,244,544),(136,216,544),(136,296,544),(136,152,544)]
    for i,shape in enumerate(unique):
        assert all(counts[shape] == 3*int(k==i) for k,counts in enumerate(shapes))
        assert rows[i][2] > 0
    print("Five metric triangle weights for radial stationarity", *[(r[2]/3).str(7) for r in rows])
    return {0:base+neutral, **signed}, 1-change, cbox, [r[2]/3 for r in rows], unique


if __name__ == "__main__":
    verify()

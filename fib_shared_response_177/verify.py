"""Exact certificates for X11; all inputs are bundle-relative and read-only.

SPDX-License-Identifier: Apache-2.0
The FIB input coefficients are consumed from fib_shared_response_176.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import flint
from flint import fmpq, fmpq_poly, fmpq_mpoly_ctx, nmod_poly, nmod_mat
import sympy as s


ROOT = Path(__file__).resolve().parent
INPUT = ROOT.parent / "fib_shared_response_176" / "global_preflight"
CTX = fmpq_mpoly_ctx.get(("u", "b", "e"), "lex")
FU, FB, FE = CTX.gens()
u, b, e, n, a = s.symbols("u b e n a")
CHECKS = []
INPUT_HASHES = {}


def check(name, condition):
    if not condition:
        raise AssertionError(name)
    CHECKS.append(name)


def read(name):
    path = INPUT / (name + ".json")
    raw = path.read_bytes()
    INPUT_HASHES[path.name] = hashlib.sha256(raw).hexdigest()
    data = json.loads(raw)
    check(name + ":variables", data["variables"] == ["u", "b", "e"])
    return CTX.from_dict({tuple(k): fmpq(v) for k, v in data["terms"]})


def conv(ex):
    return CTX.from_dict({k: fmpq(int(v.p), int(v.q))
                          for k, v in s.Poly(s.expand(ex), u, b, e, domain=s.QQ).terms()})


def univ(poly, variable=1):
    check("univariate", all(all(z == 0 for i, z in enumerate(k) if i != variable)
                            for k in poly.to_dict()))
    return fmpq_poly([poly.to_dict().get(tuple(j if i == variable else 0
                                             for i in range(3)), 0)
                      for j in range(int(poly.degrees()[variable]) + 1)])


def lc(poly, variable):
    degree = poly.degrees()[variable]
    return CTX.from_dict({tuple(0 if i == variable else z for i, z in enumerate(k)): v
                          for k, v in poly.to_dict().items() if k[variable] == degree})


def algebra():
    theta, xi = s.Rational(3, 7), s.Rational(-35, 29)
    h, ell = u * (u + 1), u*u + 3*u - 4
    al, be = 3*h - 10, 8*u*u + 18*u - 30
    phi, delta = al*a + be*n - 23*h, n*h + a*n*ell - 4*a*h
    Q = a*a-a*n-2*(u+1)*a-n*n+2*(2*u+1)*n-u*u+u+1
    T = (theta*n-a)*phi-5*xi*delta
    A = (10-3*b)*u*u+(8-3*b)*u+10*b-32
    B = (10-8*b)*u*u+(19-18*b)*u+30*b-36
    C = (23*b-29)*h
    ps6 = (n+2*a)*(10*u*u+19*u-36)-10*a*ell-29*h
    ps7 = 4*(2*n+3*a)*(20*u*u+34*u-70)-105*a*ell-381*h
    check("child line", s.expand(ps6-b*phi-(A*a+B*n+C)) == 0)
    qm = s.Matrix([[1, -s.Rational(1, 2), -u-1],
                   [-s.Rational(1, 2), -1, 2*u+1], [-u-1, 2*u+1, -u*u+u+1]])
    tm = s.Matrix([[-al, (theta*al-be-5*xi*ell)/2, (23+20*xi)*h/2],
                   [(theta*al-be-5*xi*ell)/2, theta*be, -(23*theta+5*xi)*h/2],
                   [(23+20*xi)*h/2, -(23*theta+5*xi)*h/2, 0]])
    vec = s.Matrix([A, B, C])
    rho = (vec.T*qm.adjugate()*vec)[0]
    sigma = (vec.T*tm.adjugate()*vec)[0]
    tau = (vec.T*((qm+tm).adjugate()-qm.adjugate()-tm.adjugate())*vec)[0]
    R = read("R")
    check("VI24 adjugate R", conv(tau*tau-4*rho*sigma) == R)
    qa = s.Poly(s.cancel(A*A*Q.subs(a, (-B*n-C)/A)), n)
    ta = s.Poly(s.cancel(A*A*T.subs(a, (-B*n-C)/A)), n)
    q2, q1, q0 = qa.all_coeffs()
    t2, t1, t0 = ta.all_coeffs()
    syl = (q2*t0-q0*t2)**2-(q2*t1-q1*t2)*(q1*t0-q0*t1)
    check("Sylvester source identity", conv(syl) == conv(A)**4*R)
    N = s.cancel(-203*(t2*q0-q2*t0)/A)
    D = s.cancel(-203*(q2*t1-t2*q1)/A)
    check("linear chart identity", s.expand(t2*qa.as_expr()-q2*ta.as_expr()
                                            -A*(D*n-N)/203) == 0)
    check("A input", conv(A) == read("A"))
    check("D input", conv(D) == read("D"))
    G = read("cancelled")
    check("cancelled factor", G == -conv(A)/3)
    aa, nn, zz = -(B*N+C*D), A*N, A*D
    ph = read("phi_reconstructed")
    eta = conv(4*(2*nn+3*aa)*(20*u*u+34*u-70)-105*aa*ell-381*h*zz)/G
    check("Phi sharp chart", ph*G == conv(al*aa+be*nn-23*h*zz))
    E = read("E")
    check("actual MM equation", E == FE*ph-eta)
    H, raw = read("H_factor_1"), read("H_raw")
    check("MM image resultant", R.resultant(E, 0) == raw)
    content, rem = divmod(raw, H)
    check("MM content identity", not rem and content.degrees()[2] == 0
          and H.degrees()[2] == 16 and R.degrees()[0] == 16)
    check("constant H leading coefficient", lc(H, 2).total_degree() == 0 and bool(lc(H, 2)))
    hd = read("H_main_discriminant")
    check("MM image discriminant", H.resultant(H.derivative(2), 2) == hd)
    exceptional = {"lc_R": univ(lc(R, 0)), "content": univ(content)}
    for name, poly in [("A", conv(A)), ("D", conv(D)), ("phi", ph)]:
        res = read("res_R_" + name)
        check("exceptional resultant " + name, R.resultant(poly, 0) == res)
        exceptional[name] = univ(res)
    fh, fx, fxi, fn, fd = conv(h), conv(20*u*u+34*u-70), conv(xi), conv(N), conv(D)
    U = -360*fxi*fh*fd*fn+(4*fn+(1-FU)*fd)*(210*ph-360*fxi*(fn*fx-35*fh*fd))
    W = fn*ph
    world_checks(Q, T, phi, delta, ps6, ps7, theta, xi, h)
    return R, H, hd, ph, fn, U, W, exceptional


def world_checks(Q, T, phi, delta, ps6, ps7, theta, xi, h):
    P = s.Matrix([[n+a-3*u-1, n+2*a-4*u-3], [1+2*u-a, 2+3*u-n-a]])
    M, J = s.Matrix([[0, 1], [1, 1]]), s.Matrix([[0, 1], [1, 0]])
    current, traces = s.eye(2), [1]
    for op in (M, M, J, M, M, M, M):
        current = op*current
        traces.append(s.expand(s.trace(current*P)))
    check("source trace and determinant", s.trace(P) == 1 and s.expand(P.det()-Q) == 0)
    check("literal MMJMMMM traces", traces == [1, u, u+1, n, a, n+a, n+2*a, 2*n+3*a])
    current = J*M**3*J*M**2
    check("literal MMJMMMJ trace", s.expand(s.trace(current*P)) == 4*n-u+1)
    d, t = s.symbols("d t")
    values = [-d, -d/u, -d/(u+1), (6-d)/n]
    gain = s.interpolate(list(enumerate(values)), t)
    g6 = 120/n+d*((10*u*u+19*u-36)/h-20/n)
    g7 = 210/n+d*((20*u*u+34*u-70)/h-35/n)
    check("whole cubic interpolation", s.degree(gain, t) <= 3 and all(
        s.cancel(gain.subs(t, i)-v) == 0 for i, v in enumerate(values)))
    check("clock six gain", s.cancel(gain.subs(t, 6)-g6) == 0)
    check("clock seven gain", s.cancel(gain.subs(t, 7)-g7) == 0)
    ep = d*delta-24*h*(theta*n-a)
    check("parent clock four", s.cancel(a*gain.subs(t, 4)+d-24*theta-ep/(n*h)) == 0)
    check("parent clock five", s.cancel((n+a)*gain.subs(t, 5)+d-60*(1+theta+xi)
                                        -5*ep/(2*n*h)-(d*phi-120*h*xi)/(2*h)) == 0)
    check("parent T recovery", s.cancel((ep.subs(d, 120*h*xi/phi)*phi)+24*h*T) == 0)
    th = a/n+d*delta/(24*h*n)
    check("actual sixth M", s.cancel(d+(n+2*a)*g6-120*(1+2*th)-d*ps6/h) == 0)
    check("actual terminal MM", s.cancel(d+(2*n+3*a)*g7-210*(2+3*th)-d*ps7/(4*h)) == 0)
    # Phi sharp=-3 D Phi and n=N/D. This identity uses independent symbols,
    # so it holds before imposing any source equation or positivity filter.
    NN, DD, PP = s.symbols("NN DD PP")
    jnum = 120*h*xi*n+(4*n-u+1)*(210*phi+120*xi*(n*(20*u*u+34*u-70)-35*h))
    check("actual terminal MJ", s.cancel((d+(4*n-u+1)*g7).subs(d, 120*h*xi/phi)-jnum/(n*phi)) == 0)
    raw_j = (jnum/(n*phi)).subs({n: NN/DD, a: (PP/(-3*DD)-(8*u*u+18*u-30)*NN/DD+23*h)/(3*h-10)})
    uj = -360*xi*h*DD*NN+(4*NN+(1-u)*DD)*(210*PP-360*xi*(NN*(20*u*u+34*u-70)-35*h*DD))
    check("actual MJ chart U/W", s.cancel(raw_j-uj/(NN*PP)) == 0)


def valuation(poly, factor):
    check("nonzero resultant for valuation", bool(poly))
    result = 0
    while not poly % factor:
        poly = poly // factor
        result += 1
    return result


def sturm(poly):
    chain = [poly, poly.derivative()]
    while chain[-1]:
        rem = -(chain[-2] % chain[-1])
        if not rem:
            break
        chain.append(rem / abs(rem.leading_coefficient()))
    positive = [1 if p.leading_coefficient() > 0 else -1 for p in chain]
    negative = [v*(-1)**p.degree() for v, p in zip(positive, chain)]
    variations = lambda signs: sum(x != y for x, y in zip(signs, signs[1:]))
    return variations(negative)-variations(positive)


def rational_specialization(poly, bv):
    out = [fmpq(0)] * (int(poly.degrees()[0])+1)
    for k, value in poly.to_dict().items():
        out[k[0]] += value*bv**k[1]
    return fmpq_poly(out)


def modular(poly, bv, variable=0):
    out = [0] * (int(poly.degrees()[variable])+1)
    for k, value in poly.to_dict().items():
        coefficient = int(value.numerator)*pow(int(value.denominator), -1, 101) % 101
        out[k[variable]] = (out[k[variable]] + coefficient*(pow(bv, k[1], 101)
                                                         if variable == 0 else 1)) % 101
    return nmod_poly(out, 101)


def modular_rational(poly):
    return nmod_poly([int(v.numerator)*pow(int(v.denominator), -1, 101) % 101
                      for v in poly], 101)


def coeffs(poly):
    return [int(v) for v in poly]


def det_mod(rows):
    """Independent integer Gaussian evaluation of the concrete cyclic matrix."""
    mat = [row[:] for row in rows]
    det = 1
    for j in range(len(mat)):
        pivot = next((i for i in range(j, len(mat)) if mat[i][j] % 101), None)
        if pivot is None:
            return 0
        if pivot != j:
            mat[j], mat[pivot] = mat[pivot], mat[j]
            det = -det
        value = mat[j][j] % 101
        det = det*value % 101
        inv = pow(value, -1, 101)
        for i in range(j+1, len(mat)):
            scale = mat[i][j]*inv % 101
            mat[i] = [(x-scale*y) % 101 for x, y in zip(mat[i], mat[j])]
    return det % 101


def separation(label, r, numerator, denominator):
    check(label+": degree", r.degree() in (14, 16))
    r = r/r.leading_coefficient()
    gd, ds, dt = denominator.xgcd(r)
    gs, ss, st = r.xgcd(r.derivative())
    check(label+": denominator Bezout", gd == 1 and ds*denominator+dt*r == 1)
    check(label+": squarefree Bezout", gs == 1 and ss*r+st*r.derivative() == 1)
    q = numerator*ds % r
    m, columns, power = r.degree(), [], nmod_poly([1], 101)
    for j in range(m):
        columns.append([int(power[i]) for i in range(m)])
        power = power*q % r
    rows = [[columns[j][i] for j in range(m)] for i in range(m)]
    determinant = det_mod(rows)
    check(label+": cyclic determinant", determinant != 0 and determinant == int(nmod_mat(rows, 101).det()))
    return {"prime": 101, "source_monic": coeffs(r), "U": coeffs(numerator),
            "W": coeffs(denominator), "W_bezout": [coeffs(ds), coeffs(dt)],
            "squarefree_bezout": [coeffs(ss), coeffs(st)], "q": coeffs(q),
            "cyclic_columns": columns, "cyclic_determinant": determinant}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="write the mathematical certificate beside this program")
    args = parser.parse_args()
    check("Python >=3.12", sys.version_info >= (3, 12))
    check("pinned arithmetic", s.__version__ == "1.14.0" and flint.__version__ == "0.9.0")
    R, H, hd, ph, N, U, W, exceptional = algebra()
    factors = [univ(read("support_factor_"+str(i))) for i in range(9)]
    for i, f in enumerate(factors):
        unit, fac = f.factor()
        check("irreducible F"+str(i), len(fac) == 1 and fac[0][1] == 1 and fac[0][0].degree() == f.degree())
    support = univ(read("collision_support"))
    product = fmpq_poly([1])
    for f in factors:
        product *= f
    check("complete MM support", support == product/product.leading_coefficient()
          and support.gcd(support.derivative()).degree() == 0)
    counts = [sturm(f) for f in factors]
    check("support real counts", counts == [1, 2, 2, 2, 4, 20, 24, 8, 49])
    sr = univ(R.resultant(R.derivative(0), 0))
    ihr = univ(hd)
    f44 = factors[7]
    units44 = {name: poly.gcd(f44).degree() for name, poly in exceptional.items()}
    check("F44 all chart units", all(degree == 0 for degree in units44.values()))
    check("F44 simple source and image discriminants", valuation(sr, f44) == valuation(ihr, f44) == 1)
    units_regular = {str(i): {name: poly.gcd(factors[i]).degree() for name, poly in exceptional.items()}
                     for i in (0, 8)}
    check("F0 F91 chart units", all(v == 0 for row in units_regular.values() for v in row.values()))
    check("F0 F91 source reduced", all(sr.gcd(factors[i]).degree() == 0 for i in (0, 8)))
    r0 = rational_specialization(R, fmpq(29, 23))
    w0 = rational_specialization(W, fmpq(29, 23))
    n0 = rational_specialization(N, fmpq(29, 23))
    zero_source = fmpq_poly([-1, -1, 1])
    check("F0 illegal zero n exactly", r0.gcd(n0) == zero_source and r0.gcd(w0) == zero_source)
    check("F0 Phi has no zero", r0.gcd(rational_specialization(ph, fmpq(29, 23))).degree() == 0)
    p0 = r0//zero_source
    check("F0 remaining reduced degree14", p0.degree() == 14 and p0.gcd(p0.derivative()).degree() == 0
          and p0.gcd(w0).degree() == 0)
    cert0 = separation("F0", modular_rational(p0), modular_rational(rational_specialization(U, fmpq(29, 23))),
                       modular_rational(w0))
    f91mod = modular(read("support_factor_8"), 0, variable=1)
    check("F91 primitive integral model", all(v.denominator == 1 for v in factors[8])
          and math.gcd(*(int(v.numerator) for v in factors[8])) == 1 and f91mod.degree() == 91)
    check("F91 exact residue homomorphism", f91mod(41) == 0)
    r91 = modular(R, 41)
    check("F91 preserved degree", r91.degree() == 16)
    cert91 = separation("F91", r91, modular(U, 41), modular(W, 41))
    check("advertised modular determinants", cert0["cyclic_determinant"] == 43 and cert91["cyclic_determinant"] == 47)
    # Existing X9 is a consumed result, not a new exclusion. Count its overlap
    # by exact Sturm at its rational endpoints, independently of isolation widths.
    lo, hi = fmpq(16423371109, 10**10), fmpq(1642337111, 10**9)
    chain = [factors[8], factors[8].derivative()]
    while True:
        rem = -(chain[-2] % chain[-1])
        if not rem:
            break
        chain.append(rem/abs(rem.leading_coefficient()))
    def variation_at(point):
        values = [p(point) for p in chain]
        signs = [1 if v > 0 else -1 for v in values if v]
        return sum(x != y for x, y in zip(signs, signs[1:]))
    check("X9 overlap exactly one F91 root", factors[8](lo) != 0 and factors[8](hi) != 0
          and variation_at(lo)-variation_at(hi) == 1)
    result = {"input_commit": "7c45a21dbea5cdc0974f878a96dc64281cc5db07",
              "input_sha256": INPUT_HASHES, "identities": CHECKS,
              "F44": {"degree": 44, "real_roots": 8, "source_resultant_degree": sr.degree(),
                      "image_resultant_degree": ihr.degree(), "source_valuation": 1, "image_valuation": 1,
                      "unit_gcd_degrees": units44},
              "regular_chart_unit_gcd_degrees": units_regular,
              "F0": {"b": "29/23", "b_mod_101": 32, "removed_zero_n": ["-1", "-1", "1"],
                     "remaining_degree": 14, "separation": cert0},
              "F91": {"degree": 91, "real_roots": 49, "b_mod_101": 41,
                      "source_lc_mod_101": int(r91.leading_coefficient()), "separation": cert91},
              "MM_support": {"degree": 158, "real_roots": 104},
              "joint_remaining_necessary_support": {"degree": 66, "real_roots": 54,
                  "factor_indices": [1, 2, 3, 4, 5, 6], "X9_overlap_in_F91": 1, "new_exclusions": 57},
              "arithmetic": "exact rational identities, characteristic-zero factorization/Sturm and deterministic finite-field nonzero certificates",
              "not_established": ["54 actual joint collisions", "fixed-parent joint emptiness", "whole-parent or original-root decision", "Lean verification", "independent review", "remote CI"]}
    path = ROOT / "certificate.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2)+"\n")
    else:
        check("preserved certificate matches", json.loads(path.read_text()) == result)
    print(json.dumps({"checks": len(CHECKS), "F44_MM_excluded": 8, "F0_F91_MJ_excluded": 50,
                      "MM_necessary_roots": 104, "joint_necessary_roots": 54, "new_joint_exclusions": 57}))


if __name__ == "__main__":
    main()

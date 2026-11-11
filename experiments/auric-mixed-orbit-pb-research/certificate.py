"""Exact identities and finite diagnostics for the mixed-orbit raw PB result.

Python >=3.12; dependencies numpy and sympy. Run from any directory.
The analytic universal proof is in README.md. Sampling is not that proof.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path
import sys

assert sys.version_info >= (3, 12)
import numpy as np
import sympy as sp


def exact_checks():
    C = sp.diag(1, 2, 3, 0)
    xi = sp.Matrix([-1 + 2 * sp.I, 1 - sp.I, 0, 0])
    t = sp.zeros(4)
    t[0, 3] = t[1, 3] = t[3, 0] = t[3, 1] = 1
    g = sp.I * (C * t - t * C)
    k = (C * t + t * C) / 2
    L = -(t * g + g * t)
    pure = {}
    for sigma in (1, -1):
        budgets = []
        for delta in (1, -1):
            R, S = -delta * k, t + delta * g
            ra, sa = R * xi, S * xi
            gram = (ra.H * ra)[0] * (sa.H * sa)[0] - sp.re((ra.H * sa)[0]) ** 2
            budgets.append(sp.simplify(gram))
        assert budgets == [0, 1]
        assert sp.expand((xi.H * L * xi)[0]) == 2
        pure[str(sigma)] = {"mass": 7, "left": 2, "budget_squares": [0, 1]}

    # A seven-moment separating inequality for the same fixed complete tuple.
    Rplus, Splus = -k, t + g
    assert sp.trace(Splus**2) == 14
    assert sp.trace(Rplus**2) == sp.Rational(5, 2)
    assert (Splus*xi).applyfunc(sp.expand) == sp.zeros(4, 1)
    assert (Rplus*xi).applyfunc(sp.expand) == sp.Matrix([0, 0, 0, -sp.Rational(1, 2)])

    weight = sp.diag(2, 2, sp.Rational(1, 8))
    coeff = sp.Matrix([2, 9, 3])
    slack = sp.Rational(229, 2) * weight - coeff * coeff.T
    minors = []
    for n in (1, 2, 3):
        for inds in itertools.combinations(range(3), n):
            minor = slack.extract(inds, inds).det()
            assert minor >= 0
            minors.append(str(minor))
    eps = sp.Rational(1, 100)
    m = sp.Rational(1, 4) - eps
    spectral_slack = sp.factor(m ** 2 - 229 * eps ** 2)
    assert spectral_slack > 0

    # An explicit 16-element orbit spanning family, with a symbolic determinant.
    e = sp.Symbol("epsilon", real=True, nonzero=True)
    identity = sp.eye(4) / 4
    seeds = []
    for i, j, sign in [(0, 1, 1), (0, 1, -1), (0, 2, 1), (0, 3, 1)]:
        A = sp.zeros(4)
        A[i, i], A[j, j] = sign, -sign
        seeds.append(identity + e * A)
    for i, j in itertools.combinations(range(4), 2):
        for value in (1, sp.I):
            A = sp.zeros(4)
            A[i, j], A[j, i] = value, sp.conjugate(value)
            seeds.append(identity + e * A)
    columns = []
    for s in seeds:
        A = sp.simplify((s - identity) / e)
        assert sp.expand((A ** 3 - A)) == sp.zeros(4)
        assert sp.trace(A ** 2) == 2 and sp.trace(A) == 0
        columns.append(sp.Matrix([s[i, i] for i in range(4)] +
                                 [f(s[i, j]) for i, j in itertools.combinations(range(4), 2)
                                  for f in (sp.re, sp.im)]))
    det = sp.factor(sp.Matrix.hstack(*columns).det())
    assert det != 0

    # Full-product symbolic identities for every supported Hermitian x and
    # both complex star legs. The six-square check below uses these pairings.
    xd = sp.symbols("x0:3", real=True)
    xr01, xi01, xr02, xi02, xr12, xi12 = sp.symbols(
        "xr01 xi01 xr02 xi02 xr12 xi12", real=True)
    qa, qb, qc, qd = sp.symbols("qa qb qc qd", real=True)
    X = sp.Matrix([[xd[0], xr01+sp.I*xi01, xr02+sp.I*xi02, 0],
                   [xr01-sp.I*xi01, xd[1], xr12+sp.I*xi12, 0],
                   [xr02-sp.I*xi02, xr12-sp.I*xi12, xd[2], 0],
                   [0, 0, 0, 0]])
    ts = sp.zeros(4)
    ts[0,3], ts[1,3] = qa+sp.I*qb, qc+sp.I*qd
    ts[3,0], ts[3,1] = qa-sp.I*qb, qc-sp.I*qd
    rs, ks = (C*X+X*C)/2, (C*ts+ts*C)/2
    lambda_s = sp.Symbol("lambda_s", real=True)
    for sigma in (1, -1):
        hs = sigma*sp.I*(C*X-X*C)/2
        js = sigma*sp.I*(C*ts-ts*C)/2
        for aa, bb in ((X,hs),(hs,rs),(ts,js),(ks,js),
                       (X,ts),(X,ks),(X,js),(hs,ts),(hs,ks),
                       (hs,js),(rs,ts),(rs,ks),(rs,js)):
            assert sp.expand(sp.trace(aa*bb)) == 0
        assert sp.expand(sp.trace(ks**2-js**2)) == 0
        Ls = -(X*hs+hs*X)+(X*ks+ks*X)-(ts*rs+rs*ts) \
             - lambda_s*(ts*js+js*ts)
        assert sp.expand(sp.trace(Ls)) == 0

    # Abstract real Gram/exterior decomposition, using only the stated
    # orthogonalities. This is a symbolic polynomial identity, not sampling.
    xx, hh, rr, tt, kk, jj, xr, tk, lam = sp.symbols("xx hh rr tt kk jj xr tk lam", real=True)
    H = sp.Matrix([[xx, 0, xr, 0, 0, 0],
                   [0, hh, 0, 0, 0, 0],
                   [xr, 0, rr, 0, 0, 0],
                   [0, 0, 0, tt, tk, 0],
                   [0, 0, 0, tk, kk, 0],
                   [0, 0, 0, 0, 0, jj]])
    total = 0
    for delta in (1, -1):
        a = sp.Matrix([1, delta, 0, 0, -delta, 0])
        b = sp.Matrix([0, 0, delta, 1, 0, delta * lam])
        total += (a.T * H * a)[0] * (b.T * H * b)[0] - (a.T * H * b)[0] ** 2
    even_tensor_sq = xx * tt + lam ** 2 * hh * jj + rr * kk + 2 * xr * tk
    odd_tensor_sq = hh * tt + lam ** 2 * xx * jj
    exterior = 2 * (xx * rr - xr ** 2 + hh * rr + kk * tt - tk ** 2 +
                    lam ** 2 * kk * jj + even_tensor_sq + odd_tensor_sq)
    assert sp.expand(total - exterior) == 0
    zero_coeff = sp.Matrix([2, 9])
    zero_slack = sp.Rational(85, 2) * sp.diag(2, 2) - zero_coeff * zero_coeff.T
    assert zero_slack[0, 0] == 81 and zero_slack[1, 1] == 4
    assert zero_slack.det() == 0
    return {"pure_counterwitness_by_sigma": pure,
            "seven_moment_separator": {"trace_Splus_squared": 14,
                "trace_Rplus_squared": "5/2", "pure_Splus_norm_squared": 0,
                "pure_Rplus_norm_squared": "1/4"},
            "uniform_scalar_family_cauchy_weights": ["2", "2", "1/8"],
            "uniform_scalar_family_cauchy_coefficients": [2, 9, 3],
            "cauchy_psd_all_principal_minors": minors,
            "spectral_squared_margin": str(spectral_slack),
            "symbolic_orbit_spanning_determinant": str(det),
            "full_product_symbolic_actor_star_pairings_both_sigmas": True,
            "full_product_symbolic_trace_L_all_lambda": 0,
            "abstract_exterior_polynomial_residual": 0,
            "zero_response_cauchy_psd_principal_minors": [81, 4, 0],
            "zero_response_spectral_squared_margin": str(m ** 2 - 85 * eps ** 2)}


def random_unitary(rng):
    z = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    q, r = np.linalg.qr(z)
    phase = np.diag(r)
    return q @ np.diag(phase / np.abs(phase)).conj()


def real_coordinates(a):
    return np.concatenate((a.real.ravel(), a.imag.ravel()))


def area(a, b):
    av, bv = real_coordinates(a), real_coordinates(b)
    w = np.outer(av, bv) - np.outer(bv, av)
    return float(np.linalg.norm(np.triu(w, 1)))


def hs(a):
    return float(np.linalg.norm(a))


def circ(a, b):
    return (a @ b + b @ a) / 2


def finite_diagnostics(count=1200):
    rng = np.random.default_rng(17245)
    c0 = np.diag([1., 2., 3., 0.])
    ratios, gram_margin, bound_ratio = [], [], []
    max_trace_error = 0.
    max_norm_formula_error = 0.
    max_mass_error = 0.
    masses = [1 / 3, 1., 7 / 2, 11.]
    for index in range(count):
        U, V = random_unitary(rng), random_unitary(rng)
        z = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
        x0 = np.zeros((4, 4), complex)
        x0[:3, :3] = (z + z.conj().T) / 2
        if index % 11 == 0:
            x0[:] = 0
        elif index % 11 in (1, 2, 3):
            x0[:] = 0
            x0[index % 11 - 1, index % 11 - 1] = 1
        q = rng.normal(size=2) + 1j * rng.normal(size=2)
        q *= 10. ** rng.uniform(-6, 6)
        x0 *= 10. ** rng.uniform(-6, 6)
        t0 = np.zeros((4, 4), complex)
        t0[:2, 3], t0[3, :2] = q, q.conj()
        common_scale = max(hs(x0), hs(t0))
        x0, t0 = x0 / common_scale, t0 / common_scale
        c = U @ c0 @ U.conj().T
        x, t = U @ x0 @ U.conj().T, U @ t0 @ U.conj().T
        epsilon = [1 / 100, 1 / 1000, 1 / 1000000][index % 3]
        spectrum = np.array([.25 + epsilon, .25 - epsilon, .25, .25])
        mass = masses[index % len(masses)]
        W = V @ np.diag(np.sqrt(mass * spectrum))
        max_mass_error = max(max_mass_error, abs(np.trace(W.conj().T @ W).real - mass))
        s = W @ W.conj().T
        for alpha, sigma in itertools.product((0, 2), (1, -1)):
            h = sigma * 1j * (c @ x - x @ c) / 2
            r, k = circ(c, x), circ(c, t)
            g = alpha * 1j * (c @ t - t @ c) / 2
            L = -2 * circ(x, h - k) - 2 * circ(t, r + g)
            sums, trace_areas = [], []
            for delta in (1, -1):
                R, S = x + delta * (h - k), t + delta * (r + g)
                a, b = R @ W, S @ W
                state_area = area(a, b)
                sums.append(state_area)
                trace_area = area(R, S)
                trace_areas.append(trace_area)
                gram_margin.append(state_area / mass - (.25 - epsilon) * trace_area)
                direct = hs(a) ** 2 * hs(b) ** 2 - float(np.vdot(a, b).real) ** 2
                max_norm_formula_error = max(max_norm_formula_error,
                    abs(state_area ** 2 - direct) / max(1., abs(direct)))
            T2 = sum(a * a for a in trace_areas)
            u, v, w = hs(x) * hs(h), hs(x) * hs(t), hs(t) ** 2
            assert T2 + 1e-11 >= 2 * u * u + 2 * v * v + (alpha ** 2 / 8) * w * w
            cauchy_bound = (85 / 2) if alpha == 0 else (229 / 2)
            assert hs(L) ** 2 <= cauchy_bound * T2 + 1e-10
            lhs, rhs = abs(np.trace(s @ L).real), sum(sums)
            assert lhs <= rhs + 1e-10 * mass
            ratios.append(lhs / rhs if rhs else 0.)
            bound_ratio.append(hs(L) ** 2 / T2 if T2 else 0.)
            max_trace_error = max(max_trace_error, abs(np.trace(L)))
    assert min(gram_margin) >= -1e-10
    return {"scope": "finite diagnostics only; universal result uses analytic inequalities",
            "seed": 17245, "control_samples": count, "sigma_instances_per_sample": 2,
            "alpha_values": [0, 2],
            "raw_pb_evaluations": len(ratios), "mass_values": masses,
            "epsilon_values": [1 / 100, 1 / 1000, 1 / 1000000],
            "max_lhs_over_budget_sum": max(ratios),
            "min_gram_monotonicity_margin": min(gram_margin),
            "max_norm_L_squared_over_T_squared": max(bound_ratio),
            "analytic_upper_bound_on_that_ratio": 229 / 2,
            "max_abs_trace_L": float(max_trace_error),
            "max_relative_direct_gram_vs_exterior_error": max_norm_formula_error,
            "max_mass_error": max_mass_error}


if __name__ == "__main__":
    evidence = {"identity": "auric-mixed-orbit-pb-research",
                "pinned_revision": "17245b37e7cf1381419fb570b94b23b421efd1fd",
                "python_version": sys.version.split()[0],
                "numpy_version": np.__version__, "sympy_version": sp.__version__,
                "lean_or_kernel_validation": False,
                "exact_checks": exact_checks(),
                "finite_diagnostics": finite_diagnostics()}
    output = Path(__file__).with_name("evidence.json")
    output.write_text(json.dumps(evidence, indent=2, sort_keys=True) + "\n")
    print(json.dumps(evidence, sort_keys=True))

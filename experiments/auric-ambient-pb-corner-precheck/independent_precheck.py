"""Independent exact restriction and raw-PB checks; no producer code is executed.

Frames/basis: main c0e3af176ddfd7c2fbff8f1c5d7075ff2baf20a9,
AURIC_FULL_PB_RAW_CALIBRATION_AND_COMMON_SOURCE_TRANSPORT_CONTINUATION,
equations (48.2), (49.1). Integer Gaussian matrices are formed before reduction.
Run with Python >=3.12, numpy, python-flint and sympy. Results are calculation
evidence, not Lean certification or native-source identification.
"""
from pathlib import Path
import hashlib
import json
import math
import sys
import time

import numpy as np
from flint import nmod_mat
import sympy as sp

assert sys.version_info >= (3, 12)
ROOT = Path(__file__).resolve().parent
PRODUCER = ROOT.parent / 'auric-ambient-pb-corner-research'
P = 1009
assert all(P % d for d in range(2, math.isqrt(P) + 1))
I4 = np.eye(4, dtype=object)
Z4 = np.zeros((4, 4), dtype=object)


def gm(a, b):
    """Gaussian integer matrix product, using two independent integer parts."""
    ar, ai = a
    br, bi = b
    return ar @ br - ai @ bi, ar @ bi + ai @ br


def adj(a):
    return a[0].T, -a[1].T


def conjugate(u, a):
    return gm(gm(u, a), adj(u))


def same(a, b):
    return np.array_equal(a[0], b[0]) and np.array_equal(a[1], b[1])


def basis():
    result = []
    for j in range(3):
        r = Z4.copy()
        r[j, j], r[3, 3] = 1, -1
        result.append((r, Z4.copy()))
    for j in range(4):
        for k in range(j + 1, 4):
            r, im = Z4.copy(), Z4.copy()
            r[j, k] = r[k, j] = 1
            im[j, k], im[k, j] = 1, -1
            result.extend([(r, Z4.copy()), (Z4.copy(), im)])
    return result


BASIS = basis()
FULL_BASIS = BASIS + [(I4, Z4)]
LEGS = [BASIS[j] for j in (7, 8, 11, 12)]


def coordinate_numerators(a):
    """Coordinates with common denominator 4 (before the frame denominator).

    Verify the actual full Hermitian matrix is reconstructed. This checks both
    the quotient coordinates and the retained identity input coordinate.
    """
    r, im = a
    assert same(a, adj(a))
    trace = sum(r[j, j] for j in range(4))
    nums = [4 * r[j, j] - trace for j in range(3)]
    for j in range(4):
        for k in range(j + 1, 4):
            nums.extend([4 * r[j, k], 4 * im[j, k]])
    nums.append(trace)
    rr, ri = Z4.copy(), Z4.copy()
    for n, (br, bi) in zip(nums, FULL_BASIS):
        rr += n * br
        ri += n * bi
    assert same((rr, ri), (4 * r, 4 * im))
    return np.array(nums, dtype=object)


def reduce_coords(a, denominator):
    nums = coordinate_numerators(a)
    inv = pow(4 * denominator, -1, P)
    return np.array([int(v) % P * inv % P for v in nums], dtype=np.int64)


def reconstruct_restriction():
    seed = 1
    denominators, frame_data, wedge_blocks, blocks = [], [], [], []
    indices = [(j, k) for j in range(16) for k in range(j + 1, 16)]
    c = (np.diag([1, 2, 3, 0]).astype(object), Z4)
    for f in range(48):
        u, denominator = (I4.copy(), Z4.copy()), 1
        for h in range(4):
            entries = []
            for _ in range(8):
                seed = (seed * 1664525 + 1013904223) & ((1 << 32) - 1)
                entries.append(seed % 5 - 2)
            re, im = np.array(entries[:4], dtype=object), np.array(entries[4:], dtype=object)
            d = sum(int(v) ** 2 for v in entries)
            assert 5 <= d <= 26
            denominators.append(d)
            outer_re = np.outer(re, re) + np.outer(im, im)
            outer_im = np.outer(im, re) - np.outer(re, im)
            hnum = (d * I4 - 2 * outer_re, -2 * outer_im)
            assert same(hnum, adj(hnum))
            assert same(gm(hnum, hnum), (d * d * I4, Z4))
            u = gm(u, hnum)
            denominator *= d
        dsq = denominator * denominator
        assert same(gm(adj(u), u), (dsq * I4, Z4))
        assert same(conjugate(adj(u), (I4, Z4)), (dsq * I4, Z4))
        aa = reduce_coords(conjugate(u, c), dsq)
        assert aa[-1] == 6 * pow(4, -1, P) % P
        out = np.array([reduce_coords(conjugate(adj(u), b), dsq)[3:15]
                        for b in BASIS], dtype=np.int64).T
        ww = []
        for t in LEGS:
            bb = reduce_coords(conjugate(u, t), dsq)
            assert bb[-1] == 0
            ww.append([(aa[j] * bb[k] - aa[k] * bb[j]) % P for j, k in indices])
        ww = np.array(ww, dtype=np.int64)
        wedge_blocks.append(ww)
        # Tensor contraction represents F(c,t), then the output offdiagonals.
        blocks.append(np.einsum('wi,lk->wlik', ww, out).reshape(48, 1800) % P)
        frame_data.append((u, denominator, ww))
    return (np.vstack(blocks), np.vstack(wedge_blocks),
            np.array(denominators, dtype=np.int64), frame_data)


def exact_numpy_determinant(array):
    """Modular Gaussian elimination without FLINT, using bounded int64 products."""
    a = array.copy().astype(np.int64)
    n = a.shape[0]
    assert a.shape == (n, n)
    determinant = 1
    for col in range(n):
        candidates = np.flatnonzero(a[col:, col])
        if not len(candidates):
            return 0
        row = col + int(candidates[0])
        if row != col:
            a[[row, col]] = a[[col, row]]
            determinant = -determinant
        pivot = int(a[col, col])
        determinant = determinant * pivot % P
        factors = a[col + 1:, col] * pow(pivot, -1, P) % P
        # Every multiplication is <=1008^2; no floating arithmetic/overflow.
        a[col + 1:, col + 1:] -= factors[:, None] * a[col, col + 1:]
        a[col + 1:, col + 1:] %= P
        a[col + 1:, col] = 0
    return determinant % P


def verify_actual_map(matrix, frames):
    """Direct matrix evaluation of a nontrivial full-input alternating map.

    The coefficients include identity-input pairs. Conjugate the resulting
    Hermitian output matrix, rather than using the restriction tensor again.
    """
    coeff = (np.arange(1800).reshape(120, 15) ** 2 + 17) % P
    predicted = matrix @ coeff.ravel() % P
    reads = []
    for u, denominator, ww in frames:
        # i->469 is applied only after integer adjoints have been formed.
        ur = (np.asarray(u[0] % P, dtype=np.int64)
              + 469 * np.asarray(u[1] % P, dtype=np.int64)) % P
        ud = adj(u)
        ua = (np.asarray(ud[0] % P, dtype=np.int64)
              + 469 * np.asarray(ud[1] % P, dtype=np.int64)) % P
        for values in ww @ coeff % P:
            r, im = Z4.copy(), Z4.copy()
            for value, (br, bi) in zip(values, BASIS):
                r += int(value) * br
                im += int(value) * bi
            m = (np.asarray(r % P, dtype=np.int64)
                 + 469 * np.asarray(im % P, dtype=np.int64)) % P
            rotated = ((ua @ m % P) @ ur % P) * pow(denominator ** 2, -1, P) % P
            for j in range(4):
                for k in range(j + 1, 4):
                    reads.extend([(int(rotated[j, k]) + int(rotated[k, j])) * pow(2, -1, P) % P,
                                  (int(rotated[j, k]) - int(rotated[k, j])) * pow(2 * 469, -1, P) % P])
    assert np.array_equal(predicted, np.array(reads))


def verify_certificate():
    mat, wedges, denoms, frames = reconstruct_restriction()
    assert mat.shape == (2304, 1800) and wedges.shape == (192, 120)
    cert = json.loads((PRODUCER / 'nonunital-certificate.json').read_text())
    arrays = {'nonunital-matrix.u16le': mat,
              'nonunital-input-wedges.u16le': wedges,
              'nonunital-denominators.u16le': denoms}
    hashes = {}
    for name, a in arrays.items():
        payload = a.astype('<u2').tobytes()
        supplied = (PRODUCER / name).read_bytes()
        assert payload == supplied, name
        hashes[name] = hashlib.sha256(payload).hexdigest()
        assert hashes[name] == cert['arrays'][name]['sha256']
    supplied_rows = np.frombuffer((PRODUCER / 'nonunital-pivot-rows.u16le').read_bytes(), dtype='<u2').astype(int)
    assert len(supplied_rows) == 1800 and np.all(np.diff(supplied_rows) > 0)
    mm = nmod_mat(mat.tolist(), P)
    rank = mm.rank()
    reduced, transpose_rank = mm.transpose().rref()
    independent_rows = [next(j for j in range(reduced.ncols()) if int(reduced[k, j]))
                        for k in range(transpose_rank)]
    assert np.array_equal(supplied_rows, independent_rows)
    wedge_rank = nmod_mat(wedges.tolist(), P).rank()
    assert rank == transpose_rank == 1800 and wedge_rank == 120
    # A separate elimination implementation checks the claimed oriented minor.
    determinant = exact_numpy_determinant(mat[supplied_rows, :])
    assert determinant == 20
    hashes['nonunital-pivot-rows.u16le'] = hashlib.sha256(supplied_rows.astype('<u2').tobytes()).hexdigest()
    assert hashes['nonunital-pivot-rows.u16le'] == cert['arrays']['nonunital-pivot-rows.u16le']['sha256']
    verify_actual_map(mat, frames)
    return {'shape': list(mat.shape), 'rank_mod_1009': rank,
            'wedge_shape': list(wedges.shape), 'wedge_rank_mod_1009': wedge_rank,
            'independent_numpy_minor': determinant, 'rref_rows_match': True,
            'integer_gaussian_frame_and_coordinate_checks': True,
            'direct_map_evaluation_matches': True, 'all_supplied_arrays_match': True,
            'denominator_range': [int(denoms.min()), int(denoms.max())], 'hashes': hashes}


def jordan(a, b):
    return (a * b + b * a) / 2


def K(a, b, sigma):
    return sigma * sp.I * (a * b - b * a) / 2


def raw(c, t, g, x, psi, sigma):
    h, r, k = K(c, x, sigma), jordan(c, x), jordan(c, t)
    ll = -2 * jordan(x, h - k) - 2 * jordan(t, r + g)
    values, moments = [], []
    for delta in (1, -1):
        rr, ss = x + delta * (h - k), t + delta * (r + g)
        rv, sv = rr * psi, ss * psi
        ar = sp.expand((rv.conjugate().T * rv)[0])
        bs = sp.expand((sv.conjugate().T * sv)[0])
        cross = sp.expand(sp.re((rv.conjugate().T * sv)[0]))
        values.append(sp.cancel(ar * bs - cross ** 2))
        moments.append((ar, bs, cross))
    return sp.simplify((psi.conjugate().T * ll * psi)[0]), values, moments


def verify_strict_witnesses():
    s = sp.symbols('s', real=True)
    c = sp.diag(1, 2, 3, 0)
    identities = 0
    for a in (1, 2):
        expected_b = ((80*s**6 + 12*s**4 + 6*s**2 + 1)/576 if a == 1
                      else (2*s**6 + 6*s**4 + 9*s**2 + 1)/9)
        expected_l = -(44*s**2 + 1)/12 if a == 1 else -2*(5*s**2 + 1)/3
        expected_gap = (s**2*(41+962*s**2-40*s**4)/72 if a == 1
                        else 4*s**2*(1+19*s**2-2*s**4)/9)
        for sigma in (-1, 1):
            for phase in (1, sp.I):
                t, x = sp.zeros(4), sp.zeros(4)
                t[a-1, 3], t[3, a-1] = phase, sp.conjugate(phase)
                g = K(c, t, sigma) + s * sp.eye(4)
                lam = a*(a*a-4*s*s)/(4*s*(a*a+2))
                x[a-1, a-1] = lam
                psi = sp.zeros(4, 1)
                psi[a-1], psi[3] = phase*s*(1-2*sigma*sp.I/a), 1
                lv, budgets, _ = raw(c, t, g, x, psi, sigma)
                assert sp.cancel(lv-expected_l) == 0
                assert all(sp.cancel(v-expected_b) == 0 for v in budgets)
                assert sp.cancel(lv*lv - 4*budgets[0] - expected_gap) == 0
                # The large-s witness uses an unnormalized ray, avoiding sqrt2.
                psi[a-1], psi[3] = phase, 1
                lv, budgets, _ = raw(c, t, g, sp.zeros(4), psi, sigma)
                assert sp.cancel(lv + 4*s) == 0
                assert all(sp.cancel(v-sp.Rational(a**4, 4)) == 0 for v in budgets)
                identities += 2
    y = sp.symbols('y', real=True)
    assert sp.expand(41+962*y-40*y*y - (41+sp.Rational(1919,2)*y) - 40*y*(sp.Rational(1,16)-y)) == 0
    assert sp.expand(1+19*y-2*y*y - (1+17*y) - 2*y*(1-y)) == 0
    return {'universal_symbolic_full_M4_cases': identities,
            'both_sigmas_all_four_legs': True, 'small_and_large_nonzero_s': True,
            'positive_budget_and_strict_gap_formulas': True}


def hermitian_fixture(n):
    a = sp.zeros(n)
    for j in range(n):
        a[j, j] = j - 2
        for k in range(j+1, n):
            a[j, k] = sp.Rational(j+k+1, 3) + sp.I*sp.Rational(j-k, 5)
            a[k, j] = sp.conjugate(a[j, k])
    return a


def verify_ambient():
    qs = [(1+2*sp.I, -2+3*sp.I), (sp.Rational(2,3)-sp.I, 3+sp.I/2)]
    cases = 0
    for n in (4, 5, 7):
        c = sp.diag(1, 2, 3, *([0]*(n-3)))
        for q1, q2 in qs:
            t = sp.zeros(n)
            t[0,3], t[1,3] = q1, q2
            t[3,0], t[3,1] = sp.conjugate(q1), sp.conjugate(q2)
            u = sp.zeros(n,1)
            u[0], u[1] = 2/sp.conjugate(q1), -1/sp.conjugate(q2)
            e4 = sp.eye(n)[:,3]
            assert (t*u-e4).applyfunc(sp.simplify) == sp.zeros(n,1)
            assert (jordan(c,t)*u).applyfunc(sp.simplify) == sp.zeros(n,1)
            for sigma in (-1, 1):
                # This residual obeys all three actor-row constraints while
                # retaining a general Hermitian fourth-row/spectator block.
                ee = sp.diag(2,-3,4,*([0]*(n-3)))
                if n > 4:
                    ee[3:,3:] = hermitian_fixture(n-3)
                g = K(c,t,sigma) + ee
                if n > 4:
                    z = sp.eye(n)[:,4]
                    for lam in (0,1,sp.I):
                        psi = u+lam*z
                        lv, budgets, _ = raw(c,t,g,sp.zeros(n),psi,sigma)
                        assert budgets == [0,0]
                        assert sp.simplify(lv + 2*sp.re((e4.conjugate().T*g*psi)[0])) == 0
                    # The surviving nonzero fourth-row coupling is detected.
                    assert any(raw(c,t,g,sp.zeros(n),u+lam*z,sigma)[0] != 0
                               for lam in (1,sp.I))
                # Independent full-product block budget identity with a full
                # pure functional having actor/environment/spectator entries.
                g = K(c,t,sigma) + sp.diag(2,-3,4,5,*([0]*(n-4)))
                if n > 4:
                    g[4:,4:] = hermitian_fixture(n-4)
                x = sp.zeros(n)
                x[:3,:3] = hermitian_fixture(3)
                psi = sp.Matrix([sp.Rational(j+1,2) + (j%3)*sp.I for j in range(n)])
                lv, budgets, moments = raw(c,t,g,x,psi,sigma)
                lv4, budgets4, moments4 = raw(c[:4,:4],t[:4,:4],g[:4,:4],x[:4,:4],psi[:4,:],sigma)
                extra = 0 if n == 4 else (psi[4:,:].conjugate().T * g[4:,4:]**2 * psi[4:,:])[0]
                assert sp.simplify(lv-lv4) == 0
                for b,b4,m4 in zip(budgets,budgets4,moments4):
                    assert sp.simplify(b-b4-m4[0]*extra) == 0
                cases += 1
    # Structural row-lock expansion (41.5) in N=7, including both companion
    # axes. These are regression instances; universal validity is the audited
    # dimension-free scalar argument, not inference from the sample count.
    n, sigma = 7, -1
    c = sp.diag(1,2,3,0,0,0,0)
    t = sp.zeros(n)
    t[0,3],t[1,3] = 1+2*sp.I,-2+3*sp.I
    t[3,0],t[3,1] = 1-2*sp.I,-2-3*sp.I
    g = hermitian_fixture(n)
    T, alpha = sp.symbols('T alpha',real=True)
    for j in range(3):
        k = (j+1)%3
        l, m, e = sp.zeros(n), sp.zeros(n), sp.zeros(n)
        l[j,j], m[k,k], e[j,k],e[k,j] = 1,1,1,1
        companion = sigma*sp.I*(l-m)*e
        for axis in (e,companion):
            jj = sigma*sp.I*(l-m)*axis
            nn, aa, bb = 2*l, j+1, k+1
            lam, mu = sp.Rational(aa-bb,2),sp.Rational(aa+bb,2)
            kt = jordan(c,t)
            vv = aa*jordan(nn,t)-jordan(nn,kt)
            x = T*axis+alpha*nn
            hh, rr = K(c,x,sigma),jordan(c,x)
            lhs = -2*jordan(x,hh-kt)-2*jordan(t,rr+g)
            rhs = 2*T*jordan(axis,kt)-2*mu*T*jordan(t,axis)-2*jordan(t,g)-2*alpha*(lam*T*jj+vv)
            assert (lhs-rhs).applyfunc(sp.expand) == sp.zeros(n)
            assert (l*vv*(sp.eye(n)-l)-l*(c*t-t*c)*(sp.eye(n)-l)/2).applyfunc(sp.expand) == sp.zeros(n)
    return {'full_complex_q_ambient_cases':cases, 'ambient_dimensions':[4,5,7],
            'zero_budget_remaining_leakage_detected':True,
            'full_product_spectator_budget_identity':True,
            'actor_row_symbolic_T_alpha_N7_axes':6,
            'sampling_is_not_a_universal_PB_proof':True}


def main():
    started=time.monotonic()
    cert=verify_certificate()
    print('Independent Gaussian-integer restriction and exact minor: passed',flush=True)
    witnesses=verify_strict_witnesses()
    print('Universal symbolic full-M4 strict witnesses: passed',flush=True)
    ambient=verify_ambient()
    print('Exact full-product ambient regressions: passed',flush=True)
    result={'status':'passed: independent calculation; no Lean certification',
            'source_pin':'c0e3af176ddfd7c2fbff8f1c5d7075ff2baf20a9',
            'python':sys.version,'numpy':np.__version__,'sympy':sp.__version__,
            'certificate':cert,'strict_witnesses':witnesses,'ambient':ambient,
            'elapsed_seconds':round(time.monotonic()-started,3)}
    (ROOT/'check-results.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(result,sort_keys=True))


if __name__=='__main__':
    main()

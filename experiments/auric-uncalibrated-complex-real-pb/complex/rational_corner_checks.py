"""Reconstruct the normal-form certificate factors over Q(i).

Checks the rational legality and real-coordinate reduction independently of
the prime-field frame construction, then checks the retained arrays and minors.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import sympy as s
from flint import nmod_mat

import normal_form_certificate as nf

ROOT = Path(__file__).resolve().parent


def main():
    cert = json.loads((ROOT / 'normal_form_certificate.json').read_text())
    assert s.isprime(nf.P)
    basis = []
    for j in range(4):
        a = s.zeros(4)
        a[j, j] = 1
        basis.append(a)
    for i, j in nf.PAIRS:
        a = s.zeros(4)
        a[i, j] = a[j, i] = 1
        basis.append(a)
        a = s.zeros(4)
        a[i, j], a[j, i] = s.I, -s.I
        basis.append(a)

    def mod(q):
        q = s.Rational(q)
        num, den = map(int, q.as_numer_denom())
        assert den % nf.P
        return num * pow(den, -1, nf.P) % nf.P

    def coordinates(a):
        out = [a[j, j] for j in range(4)]
        for i, j in nf.PAIRS:
            out += [(a[i, j] + a[j, i]) / 2,
                    (a[i, j] - a[j, i]) / (2 * s.I)]
        return np.array([mod(s.expand(x)) for x in out], dtype=np.int64)

    rows = []
    for vi, vec in enumerate(cert['vectors']):
        if vec is None:
            U, denominator = s.eye(4), 1
        else:
            v = s.Matrix([vec[j] + s.I * vec[j + 4] for j in range(4)])
            denominator = int((v.H * v)[0])
            U = (s.eye(4) - 2 * v * v.H / denominator).applyfunc(s.expand)
        assert denominator == cert['denominators'][vi]
        assert (U.H * U).applyfunc(s.expand) == s.eye(4)
        inputs = [coordinates(U * basis[j] * U.H) for j in nf.ACTOR]
        readouts = []
        outputs = [(U.H * b * U).applyfunc(s.expand) for b in basis]
        for j in range(3):
            readouts.append([mod((a[j, 3] + a[3, j]) / 2) for a in outputs])
            readouts.append([mod((a[j, 3] - a[3, j]) / (2 * s.I)) for a in outputs])
        frame = []
        for i in range(9):
            for j in range(i + 1, 9):
                wedge = np.array([(inputs[i][a] * inputs[j][b] -
                                   inputs[i][b] * inputs[j][a]) % nf.P
                                  for a, b in nf.WEDGE], dtype=np.int64)
                frame.extend((wedge[:, None] * np.array(rd)[None, :] % nf.P).reshape(-1)
                             for rd in readouts)
        direct, _ = nf.frame_rows(vec)
        assert np.array_equal(frame, direct)
        rows.extend(frame)
    A = np.array(rows, dtype=np.uint16)
    arrays = np.load(ROOT / 'normal_form_arrays.npz')
    T = nf.model_columns().astype(np.uint16)
    assert np.array_equal(A, arrays['constraints'])
    assert np.array_equal(T, arrays['model'])
    for key, matrix in [('constraint', A), ('model', T)]:
        assert hashlib.sha256(matrix.astype('<u2').tobytes()).hexdigest() == cert[key + '_sha256_u16le']
    minor = A[np.ix_(cert['minor_row_indices'], cert['minor_column_indices'])]
    assert minor.shape == (1768, 1768)
    assert int(nmod_mat(minor.tolist(), nf.P).det()) == 551
    # A separate small nonzero minor also certifies the model's 152 columns.
    reduced, rank = nmod_mat(T.T.tolist(), nf.P).rref()
    pivots = [next(j for j in range(T.shape[0]) if reduced[i, j]) for i in range(rank)]
    assert rank == 152
    model_det = int(nmod_mat(T[pivots, :].tolist(), nf.P).det())
    assert model_det
    result = {'rational_gaussian_corners': 33, 'rational_legality': True,
              'all_row_factors_reconstructed_over_Qi': True,
              'array_hashes_match': True, 'constraint_minor_size': 1768,
              'constraint_minor_det_mod_1009': 551,
              'model_minor_rows': pivots, 'model_minor_det_mod_1009': model_det,
              'arithmetic': 'exact rational and finite field; no floating point',
              'lean_kernel_certified': False}
    (ROOT / 'rational_corner_results.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'model_minor_rows'}, indent=2))


if __name__ == '__main__':
    main()

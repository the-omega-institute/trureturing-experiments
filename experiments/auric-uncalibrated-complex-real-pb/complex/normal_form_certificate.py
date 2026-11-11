"""Exact finite coefficient certificate for the common-corner normal form.

H4 basis: four diagonal units, then Xij,Yij in lexicographic pair order,
Yij=i(Eij-Eji). Alternating input pairs precede output coordinates.
Every frame is an actual rational Gaussian Householder unitary. Reduction
is performed only in a localization with denominators prime to P.
"""
import json
import sys

assert sys.version_info >= (3, 12)
from pathlib import Path
import hashlib
import numpy as np
from flint import nmod_mat

P = 1009
II = 469
assert II * II % P == P - 1
assert all(P % k for k in range(2, 32))
INV2 = pow(2, -1, P)
PAIRS = [(i, j) for i in range(4) for j in range(i + 1, 4)]
WEDGE = [(i, j) for i in range(16) for j in range(i + 1, 16)]
BASE = []
for i in range(4):
    a = np.zeros((4, 4), dtype=np.int64)
    a[i, i] = 1
    BASE.append(a)
for i, j in PAIRS:
    a = np.zeros((4, 4), dtype=np.int64)
    a[i, j] = a[j, i] = 1
    BASE.append(a)
    a = np.zeros((4, 4), dtype=np.int64)
    a[i, j], a[j, i] = II, -II % P
    BASE.append(a)
BASE = np.array(BASE)
EYE = np.eye(4, dtype=np.int64)
ACTOR = [0, 1, 2, 4, 5, 6, 7, 10, 11]


def coords(a):
    out = [a[i, i] % P for i in range(4)]
    for i, j in PAIRS:
        out += [(a[i, j] + a[j, i]) * INV2 % P,
                (a[i, j] - a[j, i]) * pow(2 * II, -1, P) % P]
    return np.array(out, dtype=np.int64)


def mm(a, b):
    return a @ b % P


def frame_rows(vec):
    if vec is None:
        unitary = star = EYE
        denominator = 1
    else:
        re, im = np.array(vec[:4]), np.array(vec[4:])
        denominator = int(re @ re + im @ im)
        assert denominator and denominator % P
        vp, vm = (re + II * im) % P, (re - II * im) % P
        unitary = (EYE - 2 * np.outer(vp, vm) * pow(denominator, -1, P)) % P
        star = unitary  # Householder is Hermitian before reduction.
    assert np.array_equal(mm(unitary, star), EYE)
    inputs = [coords(mm(mm(unitary, BASE[j]), star)) for j in ACTOR]
    readouts = []
    for j in range(3):
        for imag in (False, True):
            row = []
            for bk in BASE:
                a = mm(mm(star, bk), unitary)
                row.append(((a[j, 3] - a[3, j]) * pow(2 * II, -1, P)
                            if imag else (a[j, 3] + a[3, j]) * INV2) % P)
            readouts.append(row)
    readouts = np.array(readouts, dtype=np.int64)
    rows = []
    for i in range(9):
        for j in range(i + 1, 9):
            a, b = inputs[i], inputs[j]
            wedge = np.array([(a[k] * b[l] - a[l] * b[k]) % P
                              for k, l in WEDGE], dtype=np.int64)
            rows.extend((wedge[:, None] * rd[None, :] % P).reshape(-1)
                        for rd in readouts)
    return rows, denominator


def model_columns():
    cols = []
    identity_coords = coords(EYE)
    for j in range(120):
        out = np.zeros((120, 16), dtype=np.int64)
        out[j] = identity_coords
        cols.append(out.reshape(-1))
    for l in range(16):
        out = []
        for i, j in WEDGE:
            out.append(((BASE[j] if i == l else 0) -
                        (BASE[i] if j == l else 0)) if (i == l or j == l)
                       else np.zeros((4, 4), dtype=np.int64))
        cols.append(np.concatenate([coords(a) for a in out]))
    for h in BASE:
        out = []
        for i, j in WEDGE:
            a, b = BASE[i], BASE[j]
            out.append(coords((mm(mm(a, h), b) - mm(mm(b, h), a)) * II * INV2 % P))
        cols.append(np.concatenate(out))
    return np.array(cols, dtype=np.int64).T


def main():
    seed = 173
    vectors = [None]
    for _ in range(32):
        vec = []
        for _ in range(8):
            seed = (1664525 * seed + 1013904223) % (2**32)
            vec.append(seed % 5 - 2)
        vectors.append(vec)
    all_rows, denominators = [], []
    for vec in vectors:
        rows, d = frame_rows(vec)
        all_rows.extend(rows)
        denominators.append(d)
    matrix = np.array(all_rows, dtype=np.int64)
    model = model_columns()
    assert matrix.shape == (33 * 216, 1920)
    assert model.shape == (1920, 152)
    assert not np.any(matrix @ model % P)
    print('constructed exact reduction', matrix.shape, flush=True)
    model_rank = nmod_mat(model.tolist(), P).rank()
    rank = nmod_mat(matrix.tolist(), P).rank()
    print('ranks', rank, model_rank, flush=True)
    assert rank == 1768 and model_rank == 152
    # Exhibit an independently recomputable square minor, not just a rank call.
    reduced_transpose, row_rank = nmod_mat(matrix.T.tolist(), P).rref()
    row_indices = []
    for i in range(row_rank):
        row_indices.append(next(j for j in range(matrix.shape[0])
                                if reduced_transpose[i, j]))
    selected = matrix[row_indices, :]
    reduced_selected, selected_rank = nmod_mat(selected.tolist(), P).rref()
    col_indices = []
    for i in range(selected_rank):
        col_indices.append(next(j for j in range(matrix.shape[1])
                                if reduced_selected[i, j]))
    assert row_rank == selected_rank == rank
    determinant = int(nmod_mat(matrix[np.ix_(row_indices, col_indices)].tolist(), P).det())
    assert determinant == 551
    print('nonzero square minor', rank, determinant, flush=True)
    result = {
        'field_prime': P, 'imaginary_unit_image': II,
        'frame_rule': 'identity plus 32 Householders, eight integer entries per vector; LCG seed 173, recurrence 1664525*s+1013904223 mod 2^32, entry s mod 5 - 2',
        'vectors': vectors, 'denominators': denominators,
        'constraint_shape': list(matrix.shape), 'constraint_rank_mod_p': rank,
        'model_shape': list(model.shape), 'model_rank_mod_p': model_rank,
        'minor_row_indices': row_indices, 'minor_column_indices': col_indices,
        'minor_determinant_mod_p': determinant,
        'constraint_times_model_zero': True,
        'constraint_sha256_u16le': hashlib.sha256(matrix.astype('<u2').tobytes()).hexdigest(),
        'model_sha256_u16le': hashlib.sha256(model.astype('<u2').tobytes()).hexdigest(),
        'real_lift': 'All rows are reductions of real rational corner constraints. Denominators divide products of 2 and listed positive d, none divisible by 1009. Rank_mod_p=1768 implies real rank at least 1768. The 152 parameter columns have real rank at least 152 and lie in the universal common-corner kernel. Thus the full real kernel equals that normal-form image; no irrational real coefficients are excluded by the argument.',
        'normal_form': 'D(a,b)=eta(a,b)I+ell(a)b-ell(b)a+i(aHb-bHa)/2; eta real alternating, ell real linear, H Hermitian',
        'scope': 'Full real alternating Hermitian-valued coefficient space, dimension 1920; no unit-zero or covariance imposed. Finite rank certificate, not a sample of PB inequalities; no Lean certification.'
    }
    path = Path(__file__).with_name('normal_form_certificate.json')
    path.write_text(json.dumps(result, indent=2) + '\n')
    np.savez_compressed(Path(__file__).with_name('normal_form_arrays.npz'),
                        constraints=matrix.astype(np.uint16), model=model.astype(np.uint16))


if __name__ == '__main__':
    main()

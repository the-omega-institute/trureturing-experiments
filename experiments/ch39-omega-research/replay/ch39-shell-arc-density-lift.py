"""First-order lifting of all unprotected vertex cancellations under positive arc density."""
from itertools import combinations
from pathlib import Path
from runpy import run_path
from flint import arb, arb_mat

root = Path(__file__).resolve().parent
exchange = run_path(str(root / 'ch39-shell-all-current-selection.py'))
model = exchange['model']


def square_block(matrix, rows, cols):
    return arb_mat([[matrix[i, j] for j in cols] for i in rows])


def spectral_basis(matrix):
    n = matrix.nrows()
    identity = arb_mat([[int(i == j) for j in range(n)] for i in range(n)])
    roots = sorted(matrix.eig(), key=lambda v: float(v.real.mid()))
    assert all(v.imag.contains(0) for v in roots)
    values = [v.real for v in roots]
    columns = []
    for i, value in enumerate(values):
        projector = identity
        for j, other in enumerate(values):
            if i != j:
                projector = projector * (matrix - other * identity) / (value - other)
        k = max(range(n), key=lambda j: float(projector[j, j].mid()))
        assert projector[k, k] > 0
        columns.append([projector[j, k] / projector[k, k].sqrt() for j in range(n)])
    return values, arb_mat(columns).transpose()


def eigenvector_derivative(k):
    diagonal, hopping = exchange['sector_layers'](k)
    size = (7, 6, 6)[k]
    static = arb_mat(size, size)
    linear = arb_mat(size, size)
    for length in diagonal:
        density = model['density'](length)
        assert 0 < density < arb(1) / 20000
        static += diagonal[length] + hopping[length]
        linear += density * (diagonal[length] / 3 - hopping[length] / 6)
    values, frame = spectral_basis(static)
    coupling = frame.transpose() * linear * frame
    motion = arb_mat([[coupling[i, j] / (values[j] - values[i]) if i != j else 0
                       for j in range(size)] for i in range(size)])
    assert all(abs(x) < arb(1) / 32 for x in motion.entries())
    return values, motion


def disjoint_minors(matrix):
    n = matrix.nrows()
    count = 0
    for size in range(1, n // 2 + 1):
        magnitudes = []
        for left in combinations(range(n), size):
            complement = [j for j in range(n) if j not in left]
            for right in combinations(complement, size):
                if left < right:  # Transposition changes only the sign of this skew minor.
                    det = square_block(matrix, left, right).det()
                    magnitudes.append(abs(det).lower())
        assert min(magnitudes) > arb(2)**-70
        count += len(magnitudes)
        print('sector dimension', n, 'order', size, 'count', len(magnitudes),
              'smallest magnitude lower bound', min(magnitudes).str(8))
    return count


def main():
    sectors = []
    count = 0
    for k in range(3):
        values, motion = eigenvector_derivative(k)
        sectors.append(values)
        count += disjoint_minors(motion)
    assert count == 336
    assert all(abs(a - b) > arb(7) / 10000 for a, b in combinations(sum(sectors, []), 2))
    shifts = [arb(1)] + [(model['height'](i) / model['height'](1))**2 for i in range(2, 17)]
    assert all(s >= 1 and s < 23 for s in shifts)
    design = arb_mat([[1 / (s + v) - 1 / (s + sector[0]) for sector in sectors
                       for v in sector[1:]] for s in shifts])
    inverse = design.inv()
    norm = max(sum(abs(x) for x in row).upper() for row in inverse.tolist())
    assert norm < arb(2)**97
    print('First 16 residue readings, inverse maximum absolute row sum', norm.str(12))


if __name__ == '__main__':
    main()

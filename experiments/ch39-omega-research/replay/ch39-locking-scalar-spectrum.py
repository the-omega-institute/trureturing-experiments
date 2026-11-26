"""Integer Jacobian and exact mass polynomials for the color-weak vacuum."""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, product

import numpy as np
from sympy import Matrix, Poly, Symbol, factor_list

PAIRS = list(combinations(range(5), 2))
P = np.diag([int(j < 3) for i, j in PAIRS])
Q = np.eye(10, dtype=int) - P
COLOR = np.diag([1, 1, 1, 0, 0])
CONTRACTIONS = []
for k in range(5):
    wedge = np.zeros((10, 5), dtype=int)
    for i in range(5):
        if i != k:
            wedge[PAIRS.index(tuple(sorted((i, k)))), i] = 1 if i < k else -1
    CONTRACTIONS.append(wedge)


def residual_column(field, row, col, adjoint_sign):
    x = np.zeros((10, 10), dtype=int)
    y = np.zeros_like(x)
    (x if field == 0 else y)[row, col] = 1
    dx_norm = adjoint_sign*x.T @ P + P @ x
    dt = sum(k.T @ dx_norm @ k for k in CONTRACTIONS)
    residuals = [dx_norm + adjoint_sign*y.T @ Q + Q @ y,
                 x @ Q + adjoint_sign*P @ y.T,
                 adjoint_sign*x.T @ Q + P @ y, np.array([dx_norm.trace()])]
    for k in CONTRACTIONS:
        residuals.append(x @ k @ (2*COLOR - 2*np.eye(5, dtype=int)) + P @ k @ dt)
    return np.concatenate([r.ravel() for r in residuals])


def weight(row, col):
    w = tuple(int(i in PAIRS[row]) - int(i in PAIRS[col]) for i in range(5))
    return max(w, tuple(-n for n in w))


def spectrum():
    labels = list(product(range(2), range(10), range(10)))
    weights = [weight(row, col) for f, row, col in labels]
    blocks = defaultdict(list)
    for i, w in enumerate(weights):
        blocks[w].append(i)
    generators = [(-2, -2, -2, 3, 3), (0, 0, 0, 1, -1), (1, -1, 0, 0, 0)]
    t = Symbol("t")
    multiplicity = Counter()
    indices = defaultdict(lambda: [Fraction(0) for _ in generators])
    factorizations = {}
    for adjoint_sign in (1, -1):
        jac = np.column_stack([residual_column(*args, adjoint_sign) for args in labels])
        gram = jac.T @ jac
        for w, rows in blocks.items():
            other = [i for i, w2 in enumerate(weights) if w2 != w]
            assert not gram[np.ix_(rows, other)].any()
            block = gram[np.ix_(rows, rows)]
            key = (block.shape, block.tobytes())
            if key not in factorizations:
                polynomial = Matrix(block.tolist()).charpoly(t).as_expr()
                factorizations[key] = [(tuple(int(v) for v in Poly(f, t).all_coeffs()), n)
                                       for f, n in factor_list(polynomial)[1]]
            for factor, count in factorizations[key]:
                multiplicity[factor] += count
                for i, g in enumerate(generators):
                    charge = sum(a*b for a, b in zip(g, w))
                    indices[factor][i] += Fraction(count*charge**2, 2*sum(x*x for x in g))
    return multiplicity, dict(indices)


def verify():
    counts, indices = spectrum()
    assert sum(n*(len(f)-1) for f, n in counts.items()) == 400
    assert counts[1, 0] == 112
    assert indices[1, 0] == [35, 33, 32]
    assert [sum((len(f)-1)*idx[i] for f, idx in indices.items()) for i in range(3)] == [120]*3
    for f in sorted(counts):
        print(f, counts[f], tuple(map(str, indices[f])))


if __name__ == "__main__":
    verify()

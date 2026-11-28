"""Modular minors for the five-center, two-step Fourier jet constraints."""

from itertools import combinations, combinations_with_replacement, product
from math import factorial, prod

import numpy as np


def echelon(matrix, prime):
    matrix = np.asarray(matrix, dtype=np.int64) % prime
    rank = 0
    pivots = []
    for col in range(matrix.shape[1]):
        nonzero = np.flatnonzero(matrix[rank:, col])
        if not len(nonzero):
            continue
        pivot = rank + nonzero[0]
        matrix[[rank, pivot]] = matrix[[pivot, rank]]
        inverse = pow(int(matrix[rank, col]), -1, prime)
        row = matrix[rank] * inverse % prime
        matrix[rank] = row
        lower = matrix[rank + 1:]
        lower[:] = (lower - lower[:, col, None] * row) % prime
        pivots.append(col)
        rank += 1
    return matrix[:rank], pivots


def modular_rank(matrix, prime):
    _, pivots = echelon(matrix, prime)
    return len(pivots)


def modular_kernel(matrix, prime):
    rows, pivots = echelon(matrix, prime)
    free = [i for i in range(rows.shape[1]) if i not in pivots]
    basis = np.zeros((rows.shape[1], len(free)), dtype=np.int64)
    basis[free] = np.eye(len(free), dtype=np.int64)
    for row, pivot in reversed(list(zip(rows, pivots))):
        basis[pivot] = -(row @ basis) % prime
    return basis


def verify_quadratic_projection(prime):
    # Wedge with x removes all polynomial multiples of x. The finite grid
    # gives a lower bound on the rank of the corresponding polynomial map.
    x = np.array(list(product(range(-1, 3), repeat=4)), dtype=np.int64)
    tangent = np.column_stack((x, -x.sum(axis=1))) % prime
    inverse_five = pow(5, -1, prime)
    mean_square = (tangent**2).sum(axis=1) * inverse_five % prime
    square = (x**2 - mean_square[:, None]) % prime
    mean_cube = (tangent**3).sum(axis=1) * inverse_five % prime
    cube = (x**3 - mean_cube[:, None]) % prime

    def wedge(vector):
        return np.column_stack([
            (x[:, i]*vector[:, j] - x[:, j]*vector[:, i]) % prime
            for i, j in combinations(range(4), 2)
        ]).ravel()

    columns = []
    for term in (cube, mean_square[:, None]*x):
        for h, i, j in product(range(4), repeat=3):
            vector = np.zeros_like(x)
            vector[:, i] = x[:, h] * term[:, j]
            columns.append(wedge(vector))
    assert modular_rank(np.array(columns).T, prime) == 96
    for i in range(4):
        for j in range(i, 4):
            columns.append(wedge(x[:, i, None] * x[:, j, None] * square))
    assert modular_rank(np.array(columns).T, prime) == 105
    print(f"Prime {prime}: projected quartic span 96, with quadratic forms 105.")


def verify_mass_preserving_direction(prime, root, modes, powers, exponential):
    # Jets use x=ik. Specialize i and the fifth root in the same finite field.
    imaginary = next(j for j in range(2, prime) if (j*j + 1) % prime == 0)
    table = np.array([exponential[v] for v in modes]).T
    phases = {s: np.array([pow(root, s*sum(v), prime) for v in modes])
              for s in range(-2, 3)}
    jets = {s: table * phase % prime for s, phase in phases.items()}
    index = {a: i for i, a in enumerate(powers)}
    splits = []
    for alpha in powers:
        terms = []
        for i, beta in enumerate(powers):
            gamma = tuple(a-b for a, b in zip(alpha, beta))
            if min(gamma) >= 0:
                terms.append((i, index[gamma]))
        splits.append(terms)

    def coefficients(j, k):
        g = np.zeros(len(modes), dtype=np.int64)
        for place, sign in ((j, 1), (k, -1)):
            v = tuple(int(i == place) for i in range(5))
            g[modes.index(v)] = sign * (1-imaginary)
            g[modes.index(tuple(-x for x in v))] = sign * (1+imaginary)
        return g

    def contrast_rows(g, orders):
        direction = {s: jets[s] @ g % prime for s in (-1, 0, 1)}
        return np.array([
            sum(direction[0][i] * jets[s][k] - direction[s][i] * jets[0][k]
                for i, k in terms) % prime
            for s in (-1, 1) for alpha, terms in zip(powers, splits) if sum(alpha) in orders
        ])

    centers = [phases[s] for s in range(-2, 3)]
    slopes = [jets[s][j] for s in (-1, 0, 1) for j in range(1, 5)]
    assert modular_rank(centers + slopes, prime) == 17
    g = coefficients(0, 1)
    matrix = np.vstack((centers, contrast_rows(g, (1, 2, 3, 4))))
    assert np.all(matrix @ g % prime == 0)
    assert modular_rank(matrix, prime) == 60
    print(f"Prime {prime}: center and slope rank 17, direction kernel dimension one.")

    # After the linear velocity vanishes, cubics force the quadratic jets
    # into the same known ratios at all three centers.
    beta = {s: ((1-imaginary)*pow(root, s, prime) - (1+imaginary)*pow(root, -s, prime)) % prime
            for s in (-1, 0, 1)}
    quadratic = [(beta[0]*jets[s][i] - beta[s]*jets[0][i]) % prime
                 for s in (-1, 1) for i, alpha in enumerate(powers) if sum(alpha) == 2]
    kernel = modular_kernel(centers + slopes + quadratic, prime)
    assert kernel.shape == (61, 24)
    assert np.all(np.array(centers + slopes + quadratic) @ kernel % prime == 0)

    family = [contrast_rows(coefficients(j, 4), (4,)) @ kernel % prime for j in range(4)]
    pairs = list(combinations_with_replacement(range(4), 2))
    lifted = np.zeros((280, 240), dtype=np.int64)
    for j, k in product(range(4), repeat=2):
        block = pairs.index(tuple(sorted((j, k))))
        lifted[70*j:70*(j+1), 24*block:24*(block+1)] = family[k]
    assert modular_rank(lifted, prime) == 240
    print(f"Prime {prime}: four-direction quadratic coefficient rank 240 of 240.")


def verify(prime, root):
    assert root != 1 and pow(root, 5, prime) == 1
    modes = [v for v in product(range(-2, 3), repeat=5) if sum(map(abs, v)) <= 2]
    powers = sorted((a for a in product(range(5), repeat=4) if sum(a) <= 4),
                    key=lambda a: (sum(a), a))
    degrees = np.array([sum(a) for a in powers])
    denominators = [pow(prod(map(factorial, a)), -1, prime) for a in powers]
    exponential = {}
    for v in modes:
        frequency = [n - v[4] for n in v[:4]]
        exponential[v] = np.array([
            prod(n**j for n, j in zip(frequency, a)) * denominator % prime
            for a, denominator in zip(powers, denominators)
        ])
    if prime % 4 == 1:
        verify_mass_preserving_direction(prime, root, modes, powers, exponential)
    charges = np.array([sum(v) for v in modes])
    assert (len(modes), int(np.sum(charges == 0))) == (61, 21)
    differences = []
    for s in (-1, 1):
        rows = [e * (pow(root, s * sum(v), prime) - 1) for v, e in exponential.items()]
        differences.append(np.array(rows).T[degrees <= 3])
    assert modular_rank(np.vstack(differences), prime) == 40

    # Given equal scalar values and gradients, the common linear cone is x
    # after an invertible output change. Terms through degree three then
    # depend linearly on the charged vector and scalar coefficients.
    charged = charges != 0
    cubic = np.vstack(differences)[:, charged]
    low = degrees[degrees <= 3]
    quadratic = cubic[np.tile(low <= 2, 2)]
    first = cubic[np.tile(low <= 1, 2)]
    assert modular_rank(quadratic, prime) == 30
    low_index = {a: i for i, a in enumerate(powers) if sum(a) <= 3}
    block_height = len(low)
    constraints = np.zeros((4 * len(cubic) + len(first), 200), dtype=np.int64)
    constraints[4 * len(cubic):, 160:] = first
    for j in range(4):
        start = j * len(cubic)
        constraints[start:start + len(cubic), j*40:(j+1)*40] = cubic
        for alpha, i in low_index.items():
            if alpha[j] == 0:
                continue
            beta = tuple(a - int(k == j) for k, a in enumerate(alpha))
            for block in range(2):
                offset = block * block_height
                constraints[start + offset + i, 160:] = -cubic[offset + low_index[beta]]
    assert modular_rank(constraints, prime) == 190

    # Each column is a pair of light-center differences, with four output coordinates.
    def column():
        return np.zeros((2, 4, len(powers)), dtype=np.int64)

    columns = []
    for q in (-2, -1, 1, 2):
        group = [v for v in modes if sum(v) == q]
        base = exponential[group[0]]
        for v in group[1:]:
            for j in range(4):
                values = column()
                for block, s in enumerate((-1, 1)):
                    values[block, j] = (pow(root, s*q, prime) - 1) * (exponential[v] - base)
                columns.append(values)
    assert len(columns) == 144

    # R = 10 i Y is the central kinetic jet at which we exhibit a full-rank minor.
    source = np.zeros((4, len(powers)), dtype=np.int64)
    for j, l in product(range(4), range(5)):
        v = tuple(int(i == j) - int(i == l) for i in range(5))
        source[j] += exponential[v] - exponential[tuple(-i for i in v)]
    index = {a: i for i, a in enumerate(powers)}
    for block in range(2):
        for beta in (a for a in powers if sum(a) <= 3):
            values = column()
            for i, alpha in enumerate(powers):
                delta = tuple(a-b for a, b in zip(alpha, beta))
                if min(delta) >= 0:
                    values[block, :, i] = -source[:, index[delta]]
            columns.append(values)
    positive_degrees = np.array(columns)[:, :, :, degrees > 0]
    matrix = positive_degrees.reshape(len(columns), -1).T
    assert matrix.shape == (552, 214)
    assert modular_rank(matrix, prime) == 214
    print(f"Prime {prime}, root {root}: scalar rank 40, equal first jets rank 190, kinetic and amplitude rank 214.")


if __name__ == "__main__":
    for p, zeta in [(31, 2), (101, 95), (181, 42)]:
        verify(p, zeta)
        verify_quadratic_projection(p)

#!/usr/bin/env python3
"""Exact finite-chain fidelity and insertion checks; no asymptotic proof is claimed."""
from functools import lru_cache
from itertools import product
import sympy as sp

Q = sp.Rational

@lru_cache(None)
def configurations(length):
    return tuple(product(range(3), repeat=length))

@lru_cache(None)
def operators(length):
    words = configurations(length)
    lookup = {word: index for index, word in enumerate(words)}
    dimension = len(words)
    singlets, triplets = [], []
    for bond in range(length - 1):
        singlet, triplet = sp.zeros(dimension), sp.zeros(dimension)
        for column, word in enumerate(words):
            left, right = word[bond:bond + 2]
            swapped = word[:bond] + (right, left) + word[bond + 2:]
            triplet[column, column] += Q(1, 2)
            triplet[lookup[swapped], column] -= Q(1, 2)
            if left + right == 2:
                for output in range(3):
                    target = word[:bond] + (output, 2 - output) + word[bond + 2:]
                    singlet[lookup[target], column] += Q((-1) ** (left + output), 3)
        singlets.append(singlet)
        triplets.append(triplet)
    zero = sp.zeros(dimension)
    return singlets, triplets, sum(singlets, zero), sum(triplets, zero)

@lru_cache(None)
def ground_basis(length):
    hamiltonian = operators(length)[2]
    return sp.Matrix.hstack(*hamiltonian.nullspace())

def cup(length, bond):
    """Unnormalized insertion C_i=sqrt(3) F_i, using rational entries only."""
    source = configurations(length - 2)
    lookup = {word: index for index, word in enumerate(configurations(length))}
    matrix = sp.zeros(3 ** length, 3 ** (length - 2))
    for column, word in enumerate(source):
        for spin in range(3):
            target = word[:bond] + (spin, 2 - spin) + word[bond:]
            matrix[lookup[target], column] = (-1) ** spin
    return matrix

def assert_zero(matrix):
    assert matrix == sp.zeros(*matrix.shape)

def psd_exact(matrix):
    """Rational symmetric Schur elimination, including zero pivots."""
    assert matrix == matrix.T
    current = matrix
    while current.rows:
        pivot = current[0, 0]
        if pivot < 0:
            return False
        if pivot == 0:
            if any(current[0, j] != 0 for j in range(current.cols)):
                return False
            current = current[1:, 1:]
        else:
            column = current[1:, 0]
            current = current[1:, 1:] - column * column.T / pivot
    return True

def assert_le(left, right):
    assert sp.simplify(right - left).is_nonnegative is True, (left, right)

maxima = {1: sp.Integer(0)}
for length in range(2, 5):
    singlets, triplets, hamiltonian, perturbation = operators(length)
    for projector in singlets + triplets:
        assert projector == projector.T
        assert_zero(projector * projector - projector)
    for bond in range(length - 2):
        assert_zero(singlets[bond] * singlets[bond + 1] * singlets[bond]
                    - singlets[bond] / 9)
    h_spectrum = hamiltonian.eigenvals()
    gap = sp.simplify(1 - Q(2, 3) * sp.cos(sp.pi / length))
    assert gap in h_spectrum
    for energy in h_spectrum:
        if energy != 0:
            assert_le(gap, energy)
    ground_dimension = h_spectrum[0]
    assert ground_dimension == {2: 8, 3: 21, 4: 55}[length]
    compressed_eigenvalues = set()
    for magnetization in range(-length, length + 1):
        indices = [i for i, word in enumerate(configurations(length))
                   if sum(word) - length == magnetization]
        h_block = hamiltonian.extract(indices, indices)
        a_block = perturbation.extract(indices, indices)
        columns = h_block.nullspace()
        if not columns:
            continue
        kernel = sp.Matrix.hstack(*columns)
        compressed = (kernel.T * kernel).inv() * kernel.T * a_block * kernel
        compressed_eigenvalues.update(compressed.eigenvals())
    maximum = {2: sp.Integer(1), 3: sp.Integer(2),
               4: (53 + sp.sqrt(457)) / 28}[length]
    assert maximum in compressed_eigenvalues
    for eigenvalue in compressed_eigenvalues:
        assert_le(0, eigenvalue)
        assert_le(eigenvalue, maximum)
    maxima[length] = maximum
    assert_le(Q(length - 1, 2), maximum)
    assert_le(maximum, length - 1)
    print(f'L={length}: D={ground_dimension}; gamma={gap}; a_L={maximum}', flush=True)
    if length == 3:
        excited = sp.Matrix.hstack(*(hamiltonian - gap * sp.eye(27)).nullspace())
        compressed = (excited.T * excited).inv() * excited.T * perturbation * excited
        assert compressed == sp.eye(3) / 2
        mu_plus, mu_minus = -Q(1, 2), maximum - Q(1, 2)
        assert (mu_plus, mu_minus) == (-Q(1, 2), Q(3, 2))
        print(f'Table 2 OBC L=3: mu_plus={mu_plus}; mu_minus={mu_minus}', flush=True)
    if length == 4:
        assert_le(Q(52855, 100000), gap)
        assert_le(gap, Q(52865, 100000))
        print('Table 2 OBC L=4: gamma rounds to 0.5286 (exact rational enclosure)', flush=True)

for left in range(1, 4):
    for right in range(1, 5 - left):
        assert_le(maxima[left + right], maxima[left] + maxima[right] + 1)
print('Ground-space bounds and subadditivity: all m,n>=1 with m+n<=4', flush=True)

for length in range(3, 6):
    reduced = length - 2
    kernel = ground_basis(reduced)
    gram = kernel.T * kernel
    gram_inverse = gram.inv()
    _, original_triplets, _, original_a = operators(reduced)
    _, _, hamiltonian, perturbation = operators(length)
    cups = [cup(length, i) for i in range(length - 1)]
    for i, insertion in enumerate(cups):
        neighbor = sp.zeros(*insertion.shape)
        if i:
            neighbor += cups[i - 1]
        if i + 1 < len(cups):
            neighbor += cups[i + 1]
        assert_zero((hamiltonian * insertion - insertion - neighbor / 3) * kernel)
        for j, other in enumerate(cups):
            distance = abs(i - j)
            alpha = 1 if distance == 0 else Q(1, 3) if distance == 1 else 0
            assert_zero(kernel.T * insertion.T * other * kernel / 3 - alpha * gram)
            block = kernel.T * insertion.T * perturbation * other * kernel / 3
            error = block - alpha * (kernel.T * original_a * kernel)
            if distance == 0:
                if i in (0, length - 2):
                    assert_zero(error - gram / 3)
                else:
                    assert_zero(error - Q(2, 3) * gram
                                + kernel.T * original_triplets[i - 1] * kernel)
                bound = Q(2, 3)
            elif distance == 1:
                bound = Q(14, 3)
            elif distance == 2:
                bound = sp.Integer(1)
            else:
                assert_zero(block)
                bound = sp.Integer(0)
            assert psd_exact(bound ** 2 * gram - error.T * gram_inverse * error)
    print(f'L={length}: all insertion Gram/action identities and local error bounds exact', flush=True)

assert Q(2, 3) + 2 * Q(14, 3) + 2 == 12
assert 12 / Q(1, 3) == 36
print('Insertion-error constants: 12*S before normalization; 36 after gamma>=1/3', flush=True)
print('Finite checks support the mechanism; they do not establish the infinite-length limits.', flush=True)
print('ALL_CHECKS_OK')

#!/usr/bin/env python3
"""Exact computational-basis checks of arXiv:2607.03340v1, Conjecture 10.

Basis order is |x0 x1 x2 x3>, with x0 the most significant bit.
Every arithmetic operation used as evidence is in an exact algebraic number field.
"""
import sympy as S

I = S.I
c_expr = S.sqrt(2 - S.sqrt(2)) / 2
K = S.QQ.algebraic_field(I, c_expr)
z, one = K.zero, K.one
i, c = K.from_sympy(I), K.from_sympy(c_expr)
r2 = 2 * one - 4 * c * c
s = (one + r2) * c
assert c*c + s*s == one
N = 4
EDGES = ((0, 1), (0, 2), (2, 3))
BITS = tuple(tuple((j >> (N - 1 - k)) & 1 for k in range(N)) for j in range(2**N))
CUT = tuple(sum(x[a] != x[b] for a, b in EDGES) for x in BITS)


def conjugate(a):
    return K.from_sympy(S.conjugate(K.to_sympy(a)))


def norm_sq(psi):
    return sum((conjugate(a) * a for a in psi), z)


def diagonal_expectation(psi, diagonal):
    return sum((conjugate(a) * a * d for a, d in zip(psi, diagonal)), z)


def phase(psi, exp_i_gamma):
    # Source C = -sum C_e, so exp(-i gamma C)|x> = exp(+i gamma cut(x))|x>.
    return [a * exp_i_gamma**cut for a, cut in zip(psi, CUT)]


def mixer(psi, cosine, sine, denominator=1):
    # Full 16-dimensional tensor matrix exp(-i beta sum X_j).
    matrix = [[S.prod([cosine if x[k] == y[k] else -i*sine for k in range(N)])
               for y in BITS] for x in BITS]
    return [sum((m*a for m, a in zip(row, psi)), z) / denominator for row in matrix]


def state(gammas, mixers):
    psi = [one / 4] * 16
    for eg, (co, si, denominator) in zip(gammas, mixers):
        psi = mixer(phase(psi, eg), co, si, denominator)
    assert norm_sq(psi) == one
    return psi


def cut_expectation(psi, edge=(0, 1)):
    return diagonal_expectation(psi, [int(x[edge[0]] != x[edge[1]]) for x in BITS])


def check_lemma7(psi):
    # Source Lemma 7: global bit flip, and all one-qubit Z expectations zero.
    assert all(psi[j] == psi[15-j] for j in range(16))
    for k in range(N):
        assert diagonal_expectation(psi, [(-1)**x[k] for x in BITS]) == z


def fmt(a):
    return str(S.simplify(K.to_sympy(a)))


# Radius-2 rooted edge cone is the entire P4: distances from {0,1} are 0,0,1,2.
# Required open-ended continuation is 0,2,3, and vertex 3 has degree one.
degrees = [sum(v in e for e in EDGES) for v in range(N)]
assert degrees == [2, 1, 2, 1]
assert len(set((0, 2, 3))) == 3 and (0, 2) != (0, 1)
assert (0, 2) in EDGES and (2, 3) in EDGES and degrees[3] == 1
print('graph=P4; edges=01,02,23; root=01; p=2; continuation=0,2,3; leaf_degree=1')

special = state((i, -one), ((c, -s, 1), (c, s, 1)))
wanted = {0: one/2, 3: -i/2, 12: -i/2, 15: one/2}
assert all(a == wanted.get(j, z) for j, a in enumerate(special))
assert cut_expectation(special) == z
check_lemma7(special)
print('angles=(pi/2,-pi,-3*pi/8,3*pi/8) in (gamma1,gamma2,beta1,beta2) order')
print('special_state=' + str({format(j, '04b'): fmt(a) for j, a in enumerate(special) if a != z}))
print('special_root_expectation=' + fmt(cut_expectation(special)))
print('special_norm_squared=' + fmt(norm_sq(special)))
print('special_other_edge_expectations=' + str([fmt(cut_expectation(special, e)) for e in EDGES[1:]]))
print('fidelity_Lemma_7_special=PASS')

# Nearby, all four angles perturbed, no floating-point evaluation:
# gamma1 += 2 atan(1/100), gamma2 += 2 atan(1/97),
# beta1 += atan(1/101), beta2 += atan(1/103).
# exp(2i atan(1/d)) = (d+i)^2/(d^2+1).
# The mixer denominator sqrt(d^2+1)^4 = (d^2+1)^2 is rational.
def phase_shift(d):
    return (d*one+i)**2 / (d*d+1)

def perturbed_mixer(co, si, d):
    return (d*co-si, d*si+co, (d*d+1)**2)

near = state((i*phase_shift(100), -phase_shift(97)),
             (perturbed_mixer(c, -s, 101), perturbed_mixer(c, s, 103)))
check_lemma7(near)
near_value = S.simplify(K.to_sympy(cut_expectation(near)))
assert near_value.is_positive is True
assert all(a != z for a in near)
print('nearby_angle_shifts=(2*atan(1/100),2*atan(1/97),atan(1/101),atan(1/103))')
print('nearby_root_expectation=' + str(near_value))
print('nearby_root_expectation_positive=PASS; all_16_amplitudes_nonzero=PASS')
print('nearby_norm_squared=' + fmt(norm_sq(near)))
print('fidelity_Lemma_7_nearby=PASS')

# Independent source Lemma 7 distant-factorization test on P4, p=1:
# endpoints 1 and 3 have graph distance 3 > 2p.
p1 = state(((3*one+4*i)/5,), ((3*one/5,4*one/5,1),))
check_lemma7(p1)
assert diagonal_expectation(p1, [(-1)**(x[1]+x[3]) for x in BITS]) == z
print('fidelity_Lemma_7_p1_distance3_Z1Z3=0; single_Z_expectations=0; bit_flip=PASS')

uniform = state((one, one), ((one,z,1),(one,z,1)))
assert all(cut_expectation(uniform,e) == one/2 for e in EDGES)
print('zero_angle_root_expectation=1/2')
print('ALL_CHECKS_OK')

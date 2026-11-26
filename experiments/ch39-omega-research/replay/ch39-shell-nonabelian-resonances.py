"""匹配尺度以下，重弱矢量或重色矢量衰变为带荷对时的三体切割实极点。"""
from collections import Counter
from itertools import combinations, product
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx

census = run_path(str(Path(__file__).with_name("ch39-shell-charged-channel-census.py")))


def generators(indices):
    """Cartesian Pauli/Gell-Mann generators as integral matrices times a scale."""
    pair_generators = []
    for i, j in combinations(indices, 2):
        real, imag = acb_mat(5, 5), acb_mat(5, 5)
        real[i, j] = real[j, i] = 1
        imag[i, j], imag[j, i] = acb(0, -1), acb(0, 1)
        pair_generators.extend([(real, arb(1)/2), (imag, arb(1)/2)])
    i, j = indices[:2]
    diag = acb_mat(5, 5)
    diag[i, i], diag[j, j] = 1, -1
    pair_generators.insert(2, (diag, arb(1)/2))
    if len(indices) == 3:
        last = acb_mat(5, 5)
        last[i, i], last[j, j], last[indices[2], indices[2]] = 1, 1, -2
        pair_generators.append((last, 1/(2*arb(3).sqrt())))
    return pair_generators


def group_assignments(indices):
    # For a fixed source A and outgoing adjoint B the contracted group tensor
    # is the action of [T_A,T_B] on the charged representation (3,2).
    basis = generators(indices)
    counts = []
    for a, norm_a in basis:
        count, weight = 0, arb(0)
        for b, norm_b in basis:
            comm = a*b-b*a
            for color, weak in product(range(3), (3, 4)):
                root = acb_mat(5, 5)
                root[color, weak] = 1
                image = comm*root-root*comm
                for out_color, out_weak in product(range(3), (3, 4)):
                    value = image[out_color, out_weak]
                    if value != 0:  # entries before normalization are Gaussian integers
                        count += 1
                        weight += abs(value*norm_a*norm_b)**2
        assert weight.contains(3)  # C_A times the charged Dynkin index per source
        counts.append(count)
    # Both sectors admit A=T1, B=T2 with the same incoming/outgoing E_03.
    x, y = basis[0], basis[1]
    comm = (x[0]*y[0]-y[0]*x[0])*x[1]*y[1]
    factor = abs(comm[0, 0]-comm[3, 3])
    assert factor == arb(1)/2
    return counts, factor


def dot(a, b):
    return a[0]*b[0]-sum(x*y for x, y in zip(a[1:], b[1:]))


def vertex(momenta, spins):
    # Yang-Mills cubic tensor with all momenta outgoing, signature (+---).
    return sum(dot(spins[i], spins[j])*(dot(momenta[i], spins[k])-dot(momenta[j], spins[k]))
               for i, j, k in ((0, 1, 2), (1, 2, 0), (2, 0, 1)))


def residue(energy, b, a1, a2, overlap, group_factor):
    """Explicit spin component of lim (u-b^2) M/g_D^2 in matching units."""
    p = (energy**2/4-b*b).sqrt()
    diff = a1-a2
    k = ((b*b-(a1+a2)**2)*(b*b-diff*diff)).sqrt()/(2*b)
    ex, ey, ez = (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)
    # Source rest frame, outgoing spectator along -z, parent along +z.
    production = ((-energy, 0, 0, 0), (energy/2, 0, 0, -p), (energy/2, 0, 0, p))
    # Internal rest frame, decay along x and both charged spins along y.
    e1 = (b*b+a1*a1-a2*a2)/(2*b)
    decay = ((-b, 0, 0, 0), (e1, k, 0, 0), (b-e1, -k, 0, 0))
    boosted_z = (p/b, 0, 0, energy/(2*b))
    terms = [vertex(production, (ez, ex, eps_prod))*vertex(decay, (eps_dec, ey, ey))
             for eps_prod, eps_dec in ((ex, ex), (ey, ey), (boosted_z, ez))]
    assert terms[1] == terms[2] == 0
    result = abs(sum(terms))*overlap*group_factor
    assert (result-4*p*k*overlap*group_factor).contains(0) and result > 0
    return result


def verify():
    ctx.prec = 384
    weak_counts, weak_factor = group_assignments((3, 4))
    color_counts, color_factor = group_assignments((0, 1, 2))
    assert weak_factor == color_factor
    charged_matrix, neutral_matrix, radius = census["mass_matrices"]()
    charged = census["modes"](charged_matrix, radius)
    neutrals = census["modes"](neutral_matrix, radius, neutral=True)
    possible, resonances = Counter(), []
    for (i, x), (j, y), z in product(enumerate(charged), enumerate(charged), neutrals):
        # A flat source selects conjugate site profiles at its adjoint vertex
        # by orthonormality. Hence the spectator and exchanged masses both equal b.
        if 2*z.mass > 1 or x.mass+y.mass > z.mass:
            continue
        assert x.mass+y.mass < z.mass and 2*z.mass < 1
        possible[i, j] += 1
        overlap = census["coupling"](x, y, z)
        if (x.character-y.character+z.character) % 3:
            assert overlap.contains(0)
            continue
        assert overlap > 0
        energy = (1+2*z.mass)/2
        value = residue(energy, z.mass, x.mass, y.mass, overlap, weak_factor)
        resonances.append((2*z.mass, i, j, z, overlap, value))
    assert possible == {(0, 0):14, (0, 1):5, (1, 0):5}
    assert Counter((i, j) for _, i, j, _, _, _ in resonances) == {(0, 0):4, (0, 1):2, (1, 0):2}
    print("Site channels with a kinematically accessible pole / nonzero residue",
          sum(possible.values()), len(resonances))
    print("Source-wise nonzero group-component assignments, weak / color", weak_counts, color_counts)
    print("All Cartesian group assignments times ordered resonant site channels, weak / color",
          sum(weak_counts)*len(resonances), sum(color_counts)*len(resonances))
    print("X indices, external neutral (character, level), onset, overlap, residue at (1+onset)/2")
    resonances.sort(key=lambda row: row[0].mid())
    for onset, i, j, z, overlap, value in resonances:
        print(i, j, (z.character, z.level), onset.str(15), overlap.str(13), value.str(13))
    onset, i, j, z, overlap, _ = resonances[0]
    assert i == j == 0 and z.character == 0 and z.level == 2
    value = residue(arb(".8"), z.mass, charged[0].mass, charged[0].mass, overlap, weak_factor)
    assert arb(".00096437369") < value < arb(".00096437370")
    print("First pole invariant u/mu^2, onset/mu, explicit absolute residue/(g_D^2 mu^2) at E=.8 mu",
          (z.mass**2).str(15), onset.str(15), value.str(15))
    return resonances


if __name__ == "__main__":
    verify()

"""等耦合十九站点模型在匹配尺度内的带荷三体通道与交换极点。"""
from collections import Counter
from itertools import product
from pathlib import Path
from runpy import run_path
from typing import NamedTuple

from flint import acb, acb_mat, arb, ctx

budget = run_path(str(Path(__file__).with_name("ch39-shell-vector-cut-budget.py")))
soft = budget["width"]["soft"]


class Mode(NamedTuple):
    character: int
    level: int  # 1-based within the rotation character, including the neutral zero mode
    mass: arb   # in matching units
    profile: tuple


def mass_matrices():
    h = acb.zeta_zero(11431).imag
    base, neutral, charged = soft["vector_pencil"](h)
    pole, orbit = soft["pole"], soft["star"]["orbit"]
    rotation = soft["star"]["charged"]["rotation"]
    soft["star"]["charged"]["geometry"]()  # exact integer edge and clock signatures
    assert rotation[pole] == pole and {rotation[j] for j in orbit} == set(orbit)
    decrement = arb(23)/25 * arb(2)**-20 / (h*h + soft["point"]["study"]["offset"])
    for j in orbit:
        base[pole, pole] -= decrement
        base[j, j] -= decrement
        base[pole, j] += decrement
        base[j, pole] += decrement
    radius = soft["geo"]["theta"](804) / soft["geo"]["theta"](12)
    return base+charged, base+neutral, radius


def rotation_basis(character):
    groups = list(soft["star"]["charged"]["orbits"])
    if character == 0:
        groups.append((soft["pole"],))
    root = acb(-1, arb(3).sqrt())/2
    columns = []
    for group in groups:
        column = [acb(0)]*19
        for k, site in enumerate(group):
            column[site] = root**((k*character) % 3)/arb(len(group)).sqrt()
        columns.append(column)
    return acb_mat(columns).transpose()


def modes(matrix, radius, neutral=False):
    result = []
    for character in range(3):
        basis = rotation_basis(character)
        block = basis.conjugate().transpose()*matrix*basis
        values, vectors = block.eig(right=True, algorithm="rump")
        order = sorted(range(len(values)), key=lambda j: values[j].real.mid())
        assert all(values[a].real < values[b].real for a, b in zip(order, order[1:]))
        site_vectors = basis*vectors
        for level, j in enumerate(order, 1):
            val = values[j]
            assert val.imag.contains(0)
            if neutral and character == 0 and level == 1:
                assert val.real.contains(0)  # the exact constant vector
                continue
            assert val.real > 0
            raw = [site_vectors[i, j] for i in range(19)]
            norm = sum(abs(z)**2 for z in raw).sqrt()
            result.append(Mode(character, level, val.real.sqrt()/radius, tuple(z/norm for z in raw)))
    return sorted(result, key=lambda mode: mode.mass.mid())


def coupling(first, second, neutral):
    # All three outgoing site profiles. The scalar source is constant.
    return arb(19).sqrt()*abs(sum(a*b.conjugate()*z for a, b, z in
                                  zip(first.profile, second.profile, neutral.profile)))


def verify():
    ctx.prec = 384
    c_matrix, n_matrix, radius = mass_matrices()
    charged, neutral = modes(c_matrix, radius), modes(n_matrix, radius, neutral=True)
    assert len(charged) == 19 and len(neutral) == 18
    assert all(a.mass < b.mass for a, b in zip(charged, charged[1:]))
    counts, nonzero, all_open = Counter(), [], []
    for (i, first), (j, second), z in product(enumerate(charged), enumerate(charged), neutral):
        threshold = first.mass + second.mass + z.mass
        if threshold > 1:
            continue
        assert threshold < 1
        all_open.append((threshold, i, j, z))
        allowed = (first.character-second.character+z.character) % 3 == 0
        factor = coupling(first, second, z)
        if not allowed:
            assert factor.contains(0)
            continue
        assert factor > 0
        counts[i, j] += 1
        nonzero.append((threshold, i, j, z, factor))
    by_mass = sorted(nonzero, key=lambda row: row[0].mid())
    assert len(all_open) == 172 and len(by_mass) == 60 and counts[0, 0] == 6
    assert sum(count for (i, j), count in counts.items() if i != j) == 54
    # At threshold choose B,Z along x and both charged spins along y.
    # Both exchanges vanish at their source vertex and the contact tensor
    # gives 2. Thus a nonzero overlap opens the amplitude above threshold.
    print("Ordered channels per complex gauge root, kinematic / nonzero", len(all_open), len(by_mass))
    print("Lowest kinematic threshold", (2*charged[0].mass+neutral[0].mass).str(15))
    print("X indices start at 0 in increasing charged mass; Z=(rotation character, within-sector level).")
    print(" X_i  X_j  open+allowed neutral multiplicity, unordered charged pair")
    for (i, j), count in sorted(counts.items()):
        if i <= j:
            print(f" {i:3}  {j:3}  {count:3}")
    for label, rows in (("first allowed", by_mass),
                        ("first beyond X0 Xbar0", [row for row in by_mass if row[1:3] != (0, 0)]),
                        ("first with two excited charges", [row for row in by_mass if min(row[1:3]) > 0])):
        if rows:
            energy, i, j, z, factor = rows[0]
            print(label, i, j, (z.character, z.level), energy.str(15), factor.str(15))
    # For exchange of a charged state the external B vertex enforces the
    # spectator index on that propagator. Its invariant t varies between
    # (m_other+b)^2 and (sqrt(s)-m_spectator)^2.
    margins = []
    for threshold, i, j, z in all_open:
        for spectator, daughter in ((i, j), (j, i)):
            delta = (charged[daughter].mass+z.mass)**2-charged[spectator].mass**2
            assert delta > 0
            margins.append((delta, spectator, daughter, z))
    margin, i, j, z = min(margins, key=lambda row: row[0].mid())
    assert margin > arb(".00071237815")
    print("Smallest t-m_X^2 over all kinematic channels", margin.str(15), i, j, (z.character, z.level))
    below = [x.mass for x in charged if 2*x.mass < 1]
    assert len(below) == 15 and all(2*x.mass > 1 for x in charged[len(below):])
    print("Heaviest pair-producible X and lightest X+Z", below[-1].str(15),
          (charged[0].mass+neutral[0].mass).str(15))
    assert below[-1] < charged[0].mass+neutral[0].mass
    # First resonance above the specified matching range. The Dalitz line
    # enters at the pair-production threshold and lies inside just above it.
    first_parent = charged[len(below)]
    emitted = next(z for z in neutral if z.level == 1 and
                   (first_parent.character-charged[0].character+z.character) % 3 == 0)
    vertex = coupling(first_parent, charged[0], emitted)
    assert first_parent.mass > charged[0].mass+emitted.mass and vertex > arb(".06435")
    print("First genuine exchange pole above matching, onset / mu and residue's site factor",
          (2*first_parent.mass).str(15), vertex.str(15))
    return charged, neutral, by_mass


if __name__ == "__main__":
    verify()

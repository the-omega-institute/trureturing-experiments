"""Exact neighborhood counts and certified volume-tail bounds for local curvature integration."""
from collections import Counter
from itertools import combinations, permutations, product

from flint import arb, ctx


def simplex_edges(site, steps):
    vertices = [site]
    for axis in steps:
        vertex = list(vertices[-1])
        vertex[axis] = (vertex[axis] + 1) % 4
        vertices.append(tuple(vertex))
    return frozenset(tuple(sorted(pair)) for pair in combinations(vertices, 2))


def verify_neighborhood():
    chains = [
        simplex_edges(site, steps)
        for site in product(range(4), repeat=4)
        for steps in permutations(range(4))
    ]
    incidence = Counter(edge for chain in chains for edge in chain)
    assert len(chains) == 6144 and len(incidence) == 3840
    assert min(incidence.values()) == 12
    reference = chains[0]
    neighbors = [chain for chain in chains if reference & chain]
    exterior = frozenset.union(*neighbors) - reference
    assert (len(reference), len(neighbors), len(exterior)) == (10, 106, 200)
    print(f"{len(neighbors)} local factors, {len(exterior)} conditioning edges")


def verify_truncation():
    ctx.prec = 160
    shape = arb(1920)
    num_cells = arb(6144)
    conversion = num_cells * (shape.lgamma() - (shape + arb(1)/2).lgamma()).exp()
    assert arb("140.22610372") < conversion < arb("140.22610373")
    # Area <= sqrt(3)*Q/4, |2pi/n-theta| <= 2pi/3, ten hinges,
    # Q <= sqrt(96 V_sigma). Use this common bound for both tail moments.
    local_constant = 10 * arb(3).sqrt()/4 * (2*arb.pi()/3) * arb(96).sqrt()
    global_constant = local_constant * num_cells.sqrt()
    assert global_constant < 7000
    tilt = arb(2)/5
    for k in (1, 3):
        tail = (-3200*tilt - (shape + arb(k)/2)*(1-tilt).log()).exp()
        assert tail < arb(10)**-129
        error = 2 * global_constant**k * tail/(1-tail)
        bound = arb(10)**(-124 if k == 1 else -117)
        assert error < bound
        print(f"k={k}, truncated normalized moment error < {bound}")
    print("Single-simplex conversion", conversion)


if __name__ == "__main__":
    verify_neighborhood()
    verify_truncation()

"""A nearest-edge isometry for two consecutive arrivals, using rational local coins."""
from collections import Counter, defaultdict
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
from runpy import run_path

local = run_path(str(Path(__file__).with_name("ch39-local-product-contact.py")))


def uniform_coin(values):
    """Householder taking the first basis vector to the equal positive column."""
    size = len(values)
    root = isqrt(size)
    assert root**2 == size
    if size == 1:
        return values
    direction = [int(i == 0) - F(1, root) for i in range(size)]
    scale = 2*sum(x*y for x, y in zip(direction, values))/sum(x*x for x in direction)
    return [y-scale*x for x, y in zip(direction, values)]


def verify():
    sites, pole = local["lifted_sites"]()
    distance = local["geometry"]["squared"]
    nearest = [i for i, site in enumerate(sites) if distance(site, pole) == 68]
    other = [j for j in range(len(sites)) if j not in nearest]
    root = sites.index(pole)
    laplacian, _ = local["lifted_laplacian"]()
    neighbors = [{j for j in range(len(sites)) if j == i or laplacian[i,j] < 0}
                 for i in range(len(sites))]
    assert len(nearest) == 9 and len(other) == 16
    assert sum(neighbors[j].isdisjoint(nearest) for j in other) == 6

    # A path index denotes a separate port mode even when spatial edges coincide.
    paths = [("source", root, min(neighbors[root] & neighbors[j]), j) for j in other]
    for j in other:
        reachable = [i for i in nearest if neighbors[j] & neighbors[i]]
        assert len(reachable) in (3, 9)
        for i in reachable:
            mid = min(neighbors[j] & neighbors[i])
            paths.extend([("transfer", j, mid, i)]*(9//len(reachable)))
    assert len(paths) == 160
    outgoing, incoming = defaultdict(list), defaultdict(list)
    for channel, (stage, start, mid, end) in enumerate(paths):
        assert mid in neighbors[start] and end in neighbors[mid]
        outgoing[stage, start].append(channel)
        incoming[stage, end].append(channel)
    assert Counter(map(len, outgoing.values())) == {16: 1, 9: 16}
    assert Counter(map(len, incoming.values())) == {1: 16, 16: 9}

    def period(source, radial):
        transported = {}
        for (stage, node), indices in outgoing.items():
            amplitude = source if stage == "source" else radial.get(node, F(0))
            pieces = uniform_coin([amplitude] + [F(0)]*(len(indices)-1))
            transported.update(zip(indices, pieces, strict=True))
        output = {}
        for (_, node), indices in incoming.items():
            gathered = uniform_coin([transported[k] for k in indices])
            assert all(x == 0 for x in gathered[1:])
            output[node] = gathered[0]
        return output

    first = period(F(1), {})
    assert all(first[i] == F(int(i in other), 4) for i in range(len(sites)))
    for early, late in ((F(1), F(0)), (F(0), F(1)), (F(3,5), F(-4,5))):
        encoded = period(late, {j: early*v for j, v in first.items()})
        assert all(encoded[i] == (early/3 if i in nearest else late/4)
                   for i in range(len(sites)))
        assert sum(x*x for x in encoded.values()) == early**2 + late**2
        reflected_mean = F(2,25)*sum(encoded.values())
        detected = sum((reflected_mean-encoded[i])/3 for i in nearest)
        assert detected == (-7*early+24*late)/25
    print("Two nearest-edge shifts per repeated period; one shift misses six sites.")
    print("16 source paths, 144 equal transfer paths; exact rational bin isometry.")


if __name__ == "__main__":
    verify()

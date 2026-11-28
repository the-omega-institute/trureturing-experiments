"""Exact checks on labeled six-source geometries and the sharp graph order."""

from collections import deque
from itertools import combinations, permutations, product
from math import isqrt


def square_free_part(value):
    for k in range(2, isqrt(value)+1):
        while value % (k*k) == 0:
            value //= k*k
    return value


def six_sources():
    points = ((3, 6, 2, 1), (2, 3, 5, 5), (2, 5, 7, 0),
              (6, 2, 0, 5), (5, 1, 7, 5), (5, 0, 5, 3))
    pairs = list(combinations(range(len(points)), 2))
    squared = [sum((x-y)**2 for x,y in zip(points[i], points[j])) for i,j in pairs]
    radicals = [square_free_part(q) for q in squared]
    assert len(set(radicals)) == 15
    # Distinct square-free radicals are linearly independent over Q. With
    # integer masses 1 and 2 each coefficient vector is therefore an exact rate.
    rates = set()
    for tail in permutations(range(1, 6)):
        if tail[0] > tail[-1]:
            continue
        cycle = (0,) + tail
        edges = [pairs.index(tuple(sorted(e))) for e in zip(cycle, cycle[1:]+cycle[:1])]
        for masses in product((1, 2), repeat=6):
            vector = [0]*15
            for j, mass in zip(edges, masses):
                vector[j] = mass*isqrt(squared[j]//radicals[j])
            rates.add(tuple(vector))
    assert len(rates) == 3840
    print("Six Euclidean sources, two masses, distinct exponential rates:", len(rates))


def spanning_chords(n, edges):
    labels = list(range(n))
    chords = []
    for i,j in edges:
        left, right = labels[i], labels[j]
        if left == right:
            chords.append((i,j))
        else:
            labels = [left if x == right else x for x in labels]
    return chords if len(set(labels)) == 1 else None


def simultaneous_odd_length(n, edges, first, second):
    # Cover recording edge-use parity of the two specified chords. Test
    # every starting vertex because the trace sums over closed walks.
    neighbors = [[] for _ in range(n)]
    for edge in edges:
        i,j = edge
        toggle = int(edge == first)+2*int(edge == second)
        neighbors[i].append((j, toggle))
        neighbors[j].append((i, toggle))
    best = 4*n
    for start in range(n):
        queue = deque([(start, 0, 0)])
        seen = {(start, 0)}
        while queue:
            here, flags, length = queue.popleft()
            if here == start and flags == 3:
                best = min(best, length)
                break
            for there, toggle in neighbors[here]:
                state = there, flags ^ toggle
                if state not in seen:
                    seen.add(state)
                    queue.append((*state, length+1))
    return best


def graph_checks():
    for n in (4, 5):
        possible = list(combinations(range(n), 2))
        connected = largest = 0
        for pattern in product((False, True), repeat=len(possible)):
            edges = [e for e, present in zip(possible, pattern) if present]
            chords = spanning_chords(n, edges)
            if chords is None:
                continue
            connected += 1
            for a,b in combinations(chords, 2):
                length = simultaneous_odd_length(n, edges, a, b)
                largest = max(largest, length)
        assert (connected, largest) == {4:(38, 4), 5:(728, 6)}[n]
        print("Vertices, connected labeled supports, worst parity-walk length:", n, connected, largest)

    for n in range(5, 17):
        # Triangles (0,1,2) and (n-3,n-2,n-1). At n=5 vertex 2 is shared.
        first = list(combinations((0, 1, 2), 2))
        last = list(combinations((n-3, n-2, n-1), 2))
        bridges = [(i, i+1) for i in range(2, n-3)]
        edges = first + bridges + last
        assert len(bridges) == n-5
        assert simultaneous_odd_length(n, edges, first[0], last[-1]) == 2*n-4
    print("Two-triangle lower-bound family verified for 5 through 16 vertices.")


if __name__ == "__main__":
    six_sources()
    graph_checks()

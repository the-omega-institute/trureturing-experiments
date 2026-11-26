"""重超荷、弱和色矢量辐射的有序生成元群迹及交汇投影。"""
from sympy import I, Rational as Q, diag, eye, kronecker_product, sqrt, zeros


def fundamental_generators(n):
    generators = []
    for i in range(n):
        for j in range(i + 1, n):
            entry = zeros(n)
            entry[i, j] = 1
            generators.extend(((entry + entry.T)/2, (entry - entry.T)/(2*I)))
    for k in range(1, n):
        diagonal = [1]*k + [-k] + [0]*(n - k - 1)
        generators.append(diag(*diagonal)/sqrt(2*k*(k + 1)))
    return generators


def projections(vector):
    y, w, c = vector
    return (115*y - 333*w + 218*c)/218, 5*y - 3*w - 2*c


def verify():
    # A charged root maps the weak doublet into the color triplet. The weak
    # generator therefore acts in the conjugate fundamental representation.
    groups = [
        [5*eye(6)/sqrt(60)],
        [kronecker_product(eye(3), -t.T) for t in fundamental_generators(2)],
        [kronecker_product(t, eye(2)) for t in fundamental_generators(3)],
    ]
    indices = (Q(5, 2), Q(3, 2), Q(1))
    symmetric_rows, commutator_rows = [], []
    for emitter_id, emitters in enumerate(groups):
        symmetric_row, commutator_row = [], []
        for source_id, sources in enumerate(groups):
            sums = []
            for source in sources:
                assert (source**2).trace() == indices[source_id]
                symmetric_weight, commutator_weight = 0, 0
                for emitted in emitters:
                    left, right = source*emitted, emitted*source
                    symmetric, difference = (left + right)/2, (left - right)/2
                    assert (symmetric.H*difference).trace() == 0
                    symmetric_weight += (symmetric.H*symmetric).trace()
                    commutator_weight += (difference.H*difference).trace()
                    # The cubic source-heavy-heavy graph contracts the Lie
                    # bracket with the generator acting on the charged pair.
                    if source_id == emitter_id and source_id != 0:
                        bracket = left - right
                        expanded = sum(((t*bracket).trace()*t for t in sources), zeros(6))
                        assert expanded/indices[source_id] == 2*difference
                    else:
                        assert difference == zeros(6)
                sums.append((symmetric_weight, commutator_weight))
            assert all(row == sums[0] for row in sums)
            symmetric_row.append(sums[0][0])
            commutator_row.append(sums[0][1])
        symmetric_rows.append(tuple(symmetric_row))
        commutator_rows.append(tuple(commutator_row))
    assert symmetric_rows == [
        (Q(25, 24), Q(5, 8), Q(5, 12)),
        (Q(15, 8), Q(3, 8), Q(3, 4)),
        (Q(10, 3), Q(2), Q(7, 12)),
    ]
    assert commutator_rows == [(0, 0, 0), (0, Q(3, 4), 0), (0, 0, Q(3, 4))]
    for name, rows in (("anticommutator / 2", symmetric_rows), ("commutator / 2", commutator_rows)):
        print(name, "by emitted Y, W, C", rows)
    symmetric = tuple(map(sum, zip(*symmetric_rows)))
    commutator = tuple(map(sum, zip(*commutator_rows)))
    assert symmetric == (Q(25, 4), 3, Q(7, 4))
    assert commutator == (0, Q(3, 4), Q(3, 4))
    assert projections(symmetric) == (Q(405, 872), Q(75, 4))
    assert projections(commutator) == (Q(-345, 872), Q(-15, 4))
    print("Summed weights and (P, J)", symmetric, projections(symmetric),
          commutator, projections(commutator))


if __name__ == "__main__":
    verify()

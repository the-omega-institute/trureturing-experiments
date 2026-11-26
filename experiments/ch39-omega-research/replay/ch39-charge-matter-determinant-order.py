"""在 11431 的同一星图背景上判定五重态与十重态的绝对匹配投影。"""
from collections import defaultdict
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx

soft = run_path(str(Path(__file__).with_name("ch39-soft-nonspoke-vector-mode.py")))
geo, pole = soft["geo"], soft["pole"]


def nonspokes(q, height):
    phases = geo["phases"](height)
    matrix = acb_mat(19, 19)
    for i, j, d, z in geo["edges"]:
        if pole in (i, j):
            continue
        link = (z*phases.get(d, 1))**q
        matrix[i, i] += 1
        matrix[j, j] += 1
        matrix[i, j] = -link
        matrix[j, i] = -link.conjugate()
    return matrix


def coefficients(base, nonspoke):
    roots = [acb(0, 2*arb.pi()*k/20).exp() for k in range(20)]
    values = [(-base-r*nonspoke).charpoly() for r in roots]
    table = {}
    for j in range(20):
        transform = acb.dft([p[j] for p in values])
        for i in range(min(18, 19-j)+1):
            table[i, j] = transform[i].real/20
    return table


def support(i, j):
    # A phase-dependent permutation must contain an oriented cycle. A cycle
    # with a spoke needs two spokes; one without a spoke needs three other edges.
    return i >= 1 and (i >= 3 or 19-j-i >= 2)


def mass_wronskian(differences, denominator):
    wronskian = defaultdict(arb)
    for (i, j), diff in differences.items():
        for (l, k), value in denominator.items():
            if k != j and (l, k) != (0, 0):
                wronskian[i+l, j+k-1] += (k-j)*diff*value
    # At e^18 each determinant has the common factor x+18 from the pole.
    # Its contribution at e^36 therefore cancels in the mass derivative.
    assert wronskian.pop((36, 0)).contains(0)
    return wronskian


def contrast(index, weight, mass2):
    h = acb.zeta_zero(index).imag
    base = soft["vector_pencil"](h)[0]+mass2*soft["identity"]
    logs = {q:(base+weight*nonspokes(q, h)).det().real.log() for q in (1, 2, 3, 4, 6)}
    return (2*(logs[2]-logs[3])/(5*arb.pi()),
            (3*logs[4]+logs[6]-4*logs[1])/(5*arb.pi()))


def verify():
    ctx.prec = 224
    height = acb.zeta_zero(11431).imag
    base = soft["vector_pencil"](height)[0]
    tables = {q:coefficients(base, nonspokes(q, height)) for q in (1, 2, 3, 4, 6)}
    for a, b in ((3, 2), (1, 4), (1, 6)):
        differences = {}
        for pair, value in tables[a].items():
            difference = value-tables[b][pair]
            if support(*pair):
                assert difference > 0, (a, b, pair)
                differences[pair] = difference
            else:
                assert difference.contains(0)
        assert len(differences) == 185
        derivative = mass_wronskian(differences, tables[b])
        assert len(derivative) == 696 and all(v > 0 for v in derivative.values())
    five, ten = contrast(11431, arb(1), arb(0))
    assert five < -arb(".0164") and ten < -arb(".0869")
    other_five = contrast(35, arb(1), arb(0))[0]
    other_ten = contrast(1, arb(1), arb(0))[1]
    assert other_five > arb(".0107") and other_ten > arb(".0142")
    print("Three sets of 185 positive determinant and 696 positive derivative coefficients.")
    print("Signed Dirac matter contributions at 11431, five and ten", five.str(14), ten.str(14))
    print("Opposite signs at actual zeros 35 (five), 1 (ten)", other_five.str(13), other_ten.str(13))


if __name__ == "__main__":
    verify()

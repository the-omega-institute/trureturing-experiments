"""在十九个原生站位的同一坐标中，由相邻实际壳本征线的重叠平方输运振幅。"""
from flint import acb, acb_mat, ctx


def spectrum(shell, q, number):
    ctx.prec = max(ctx.prec, 320)
    height = acb.zeta_zero(number).imag
    laplacian = shell['laplacian'](q, height)
    values, vectors = laplacian.eig(right=True, algorithm='rump')
    order = sorted(range(19), key=lambda j: values[j].real.mid())
    levels = [values[j].real for j in order]
    assert all(v.imag.contains(0) for v in values)
    assert all(left < right for left, right in zip(levels, levels[1:]))
    norms = [sum(abs(vectors[i, j])**2 for i in range(19)).sqrt() for j in order]
    assert all(n > 0 for n in norms)
    basis = acb_mat([[vectors[i, j]/n for j, n in zip(order, norms, strict=True)] for i in range(19)])
    gram = basis.conjugate().transpose()*basis
    assert all(gram[i, j].contains(int(i == j)) for i in range(19) for j in range(19))
    return levels, basis


def training_values(levels, training):
    amplitudes = []
    for energy in levels:
        matches = [value for node, value in training if energy.overlaps(node)]
        assert len(matches) == 1
        amplitudes.append(matches[0])
    return amplitudes


def transport(previous_basis, current_basis, amplitudes):
    assert len(amplitudes) == 19
    overlaps = current_basis.conjugate().transpose()*previous_basis
    weights = [[abs(overlaps[j, i])**2 for i in range(19)] for j in range(19)]
    assert all(sum(row).contains(1) for row in weights)
    assert all(sum(col).contains(1) for col in zip(*weights, strict=True))
    return [sum(w*a for w, a in zip(row, amplitudes, strict=True)) for row in weights]

"""允许源极点单独重相仍不能由同一物理酉矩阵在右外壳实现完整夸克中心。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb


def verify(root):
    model = run_path(str(root/'ch39-all-spaces-common-unitary-quark-right.py'))
    h, spectra = model['certify_observables'](root)
    pole = 14
    keep = [i for i in range(19) if i != pole]
    principal = acb_mat([[h[i,j] for j in keep] for i in keep])
    # 用微小的正对角项分离重根；扣去扰动算子范数上界后仍为原谱的下界。
    step = arb('1e-7')
    for j in range(18):
        principal[j,j] += (j+1)*step
    values = principal.eig(algorithm='rump')
    assert all(z.imag.contains(0) for z in values)
    lower = [arb(0)] + [max(arb(0), (x-18*step).lower()) for x in sorted(z.real.lower() for z in values)]
    thresholds = [arb(model['TRACE']), arb(model['CHARACTER']), arb(model['CHARACTER'])]
    ranks = (2,1,1)
    margins = []
    for eig, threshold, rank in zip(spectra, thresholds, ranks, strict=True):
        total = sum(a*b.lower() for a,b in zip(lower[:len(eig)], reversed(eig), strict=True))
        assert total > threshold
        margins.append(((total-threshold)/(rank*arb(model['CAPS'][0]))).lower())
        print('允许源极点重相的加权下界',total.str(14),flush=True)
    assert min(margins) > arb('.000310')
    print('输运到附加源极点相位的共同输运族的距离下界',min(margins).str(14))


if __name__ == '__main__':
    verify(Path(__file__).resolve().parent)

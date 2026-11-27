"""允许任意共同三源换基，首零点有理制备的两部门质量层级仍不相容。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def singular_values(matrix):
    values = (matrix.conjugate().transpose()*matrix).eig(algorithm="rump")
    assert all(v.is_finite() and v.imag.contains(0) and v.real > 0 for v in values)
    return sorted((v.real.sqrt() for v in values), key=lambda x: x.mid())


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    graph = lib("three-ray-chiral-boundary")
    ctx.prec = 256
    gram = lib("native-rh-ray-pvm")["gram_cholesky"](root)
    ports = [graph["soft"]["pole"],*graph["boundary_sites"]()]
    sources = acb_mat([[gram[j][ports.index(i)] if i in ports else 0 for j in (1,2,3)]
                       for i in range(19)])
    assert sum(abs(z)**2 for z in sources.entries()) < 4
    spectrum = lib("adjacent-pvm-transport")["spectrum"]
    blocks = []
    condition = arb(1)
    error = 2*max((arb("8.685")/9)**512,(9/arb("9.6"))**512)
    for q in (-4,2):
        levels,basis = spectrum(graph,q,1)
        assert 0 < levels[0] and levels[2] < arb("8.685") and levels[3] > arb("9.6")
        light = acb_mat([[basis[i,j] for j in range(3)] for i in range(19)])
        block = light.conjugate().transpose()*sources
        blocks.append(block)
        sigma = singular_values(block)[0]
        perturbation = error/sigma
        assert 0 < perturbation < arb(7)/9
        condition *= (1+perturbation)/(arb(7)/9-perturbation)
    relative = singular_values(blocks[0]*blocks[1].inv())
    ratio = relative[-1]/relative[0]
    bound = ratio*condition
    required = arb(".8573")/arb(".00000629") * arb(".0000133")/arb(".01410")
    assert ratio < arb("5.96302") and bound < arb("9.874")
    assert required > arb("128.56") and required/bound > arb("13.02")
    print("首零点相对条件数、所有目的背景上界、质量窗下界、部门附加变形至少需要的倍数",
          ratio.str(13),bound.str(13),required.str(13),(required/bound).str(13))
    def box(mid,radius):
        return arb(mid)+arb(0,arb(radius).upper())
    t12,t23,t13,phase = (box(m,r) for m,r in ((".22704",".00082"),(".04275",".00042"),
                                                (".00377",".00008"),("1.139",".023")))
    vtd = abs(acb(t12.sin()*t23.sin())-t12.cos()*t23.cos()*t13.sin()*acb(0,phase).exp())
    mixing = arb(".8573")/arb(".0000137")*arb(".01392")/arb(".00000629")*vtd.lower()*t13.sin().lower()
    assert mixing > 4329 and mixing/bound > arb("438.5")
    threshold = (required/ratio).sqrt()
    minimal_error = (threshold*arb(7)/9-1)/(threshold+1)
    assert minimal_error > arb(".46274")
    print("含CKM的相对条件数下界、附加变形条件数积下界、只用质量时相对源误差下界",
          mixing.str(13),(mixing/bound).str(13),minimal_error.str(13))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)

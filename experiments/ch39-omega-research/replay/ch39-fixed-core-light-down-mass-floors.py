"""七项见证的全部固定核心补全上，两轻下型质量的严格下界。"""
from math import prod
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb, acb_mat, arb, ctx


def destination_trivial_floor(source, target, sites, labels, fixed, images, adj):
    source_coords = adj(source)*fixed
    rhs = adj(images)*target
    shifts = {0:".00017",5:".0013",6:".0128",9:".0208",11:".00041",15:".0023",18:".082"}
    floors = [arb(0)]*19
    for j,text in shifts.items():
        shift = arb(text)
        scaled = acb_mat([[source_coords[i,k]/(labels[i]+shift) for k in range(4)] for i in range(19)])
        metric = adj(source_coords)*scaled
        data = acb_mat([[rhs[i,j]] for i in range(4)])
        trial = metric.mid().solve(data.mid()).mid()
        floors[j] = max(arb(0), (2*(adj(trial)*data)[0,0].real
                                 -(adj(trial)*metric*trial)[0,0].real-shift).lower())
    basis = acb_mat([[1,0],[0,1/arb(3).sqrt()],[0,1/arb(3).sqrt()],[0,1/arb(3).sqrt()]])
    four = acb_mat([[target[i,j] for j in range(19)] for i in sites])
    rows = adj(basis)*four
    block = acb_mat([[rows[i,j]*floors[j] for j in range(19)] for i in range(2)])*adj(rows)
    level = arb(".0035")
    block -= acb_mat([[level,0],[0,level]])
    assert block[0,0].real>0 and block.det().real>0
    return level


def free_trivial_floors(root, ranktwo, trivial, spectrum, shells, n, source, target, sites, labels):
    joint = run_path(str(root / "ch39-one-generator-up-and-gravity.py"))
    _, (a, _), gev, (low, high, _) = joint["verify"](root, return_parameters=True)
    _, six_source = spectrum(shells, 6, n)
    _, six_target = spectrum(shells, 6, n+1)
    seed = acb_mat([[acb(x,y) for x,y in zip(re,im,strict=True)]
                    for re,im in zip(ranktwo["REAL"],ranktwo["IMAG"],strict=True)])
    adj = ranktwo["adj"]
    kp = trivial * seed * (arb(1)/10000)
    projection = kp * (adj(kp)*kp).inv() * adj(kp)
    _, offset, small, big, tangent, length = joint["mode_parameters"](projection, trivial, six_source, six_target)
    e = ((gev/arb("1.22089e19")-low)/(high-low)).sqrt()
    linear = offset + small*e + big*a
    radius = (length-(adj(linear)*linear)[0,0].real).sqrt()
    center = offset + small*e.mid() + big*a.mid()
    image = center + tangent*(length-(adj(center)*center)[0,0].real).sqrt()
    for var, direction in ((e,small),(a,big)):
        image += (direction-tangent*((adj(linear)*direction)[0,0].real/radius))*(var-var.mid())
    initial = (ranktwo["eye"](19)-projection)*joint["col"](six_source,17)
    first = [joint["col"](kp,0),joint["col"](kp,1)]
    fixed = joint["side_by_side"]([*first, initial, image])
    images = joint["side_by_side"]([*(v*acb(0,1) for v in first), image, initial])
    level = destination_trivial_floor(source,target,sites,labels,fixed,images,adj)
    cols = adj(trivial)*joint["side_by_side"]([*first,initial])
    free = ranktwo["eye"](7)-cols*(adj(cols)*cols).inv()*adj(cols)
    last = free*adj(trivial)*image
    free -= last*adj(last)*(1/(adj(last)*last)[0,0].real)
    overlaps = adj(source)*trivial
    inputs = acb_mat([[overlaps[i,j]*labels[i] for j in range(7)] for i in range(19)])
    source_block = adj(overlaps)*inputs
    mid = free*source_block*free
    # 浮点本征矢只提供试探空间。每个结论由 Arb 正定检验给出。
    arr = np.array([[complex(z.mid()) for z in row] for row in mid.tolist()])
    _, vectors = np.linalg.eigh((arr+arr.conj().T)/2)
    floors = [arb(".26"),arb(".052"),arb(".00025")]
    for k, floor in enumerate(floors, 1):
        guess = acb_mat([[acb(float(z.real),float(z.imag)) for z in row] for row in vectors[:,-k:]])
        test = free*guess
        coords = overlaps*test
        scaled = acb_mat([[coords[i,j]*labels[i] for j in range(k)] for i in range(19)])
        difference = adj(coords)*scaled - (adj(test)*test)*floor
        for j in range(1,k+1):
            minor = acb_mat([[difference[i,l] for l in range(j)] for i in range(j)])
            assert minor.det().real>0, (k,j,minor.det())
    return floors, level


def character_sum_bound(root):
    ranktwo = run_path(str(root / "ch39-rank-two-generator-up-quark.py"))
    shells = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    spectrum = run_path(str(root / "ch39-adjacent-pvm-transport.py"))["spectrum"]
    train = run_path(str(root / "ch39-common-embedding-central-right-ckm.py"))["verify_nodes"]
    orbits = run_path(str(root / "ch39-up-quark-arbitrary-site-phase-floor.py"))["ORBITS"]
    n = ranktwo["N"] + 1
    _, labels, _ = train(root, ((n, ".0294"),), ".1", (), ((6,7,10,13), (6,7,10,12)))
    ctx.prec = 320
    _, source = spectrum(shells, 2, n)
    _, target = spectrum(shells, 2, n+1)
    adj = ranktwo["adj"]
    trivial = acb_mat([[1/arb(len(o)).sqrt() if i in o else 0 for o in orbits] for i in range(19)])
    overlaps = adj(trivial) * source
    costs_source = [labels[n,2][j].lower() for j in range(19)
                    if sum(abs(overlaps[i,j])**2 for i in range(7)) < arb(".01")]
    assert len(costs_source) == 12
    sites = (shells["soft"]["pole"], *shells["boundary_sites"]())
    floors,level = free_trivial_floors(root,ranktwo,trivial,spectrum,shells,n,source,target,sites,labels[n,2])
    costs_source += floors
    four = acb_mat([[target[i,j] for j in range(19)] for i in sites])
    omega = acb(-arb(1)/2, arb(3).sqrt()/2)
    characters = acb_mat([[0,1,t,t**2] for t in (omega,omega.conjugate())]) * (1/arb(3).sqrt())
    four_character = characters * four
    cost = [sum(abs(four_character[i,j])**2 for i in range(2)) for j in range(19)]
    positive = sorted(z.lower() for z in cost if z > arb(".00038"))
    assert len(positive) == 12
    # 去掉四维平凡核心后，目的代价有三个零本征值。迹重排取两列反序配对。
    return sum(v*w for v,w in zip(sorted(costs_source,reverse=True),[arb(0)]*3+positive,strict=True)),level


def bounds(reference, masses, char_sum, level):
    probs = [prod(x-y for y in masses)/prod(x-z for j,z in enumerate(reference) if j != i)
             for i,x in enumerate(reference)]
    assert all(p > 0 for p in probs)
    p_e, p_f = probs[:2]
    p_trivial = 1 - p_e - p_f
    # 固定字符总和时，最小压缩本征值为凹函数，在两端取其下界。
    def secular(z, p_zero, p_other):
        return -p_zero/z + p_other/(char_sum-z) + p_trivial/(level-z)

    lo, hi = arb(0), char_sum.lower()
    for _ in range(32):
        middle = (lo+hi)/2
        if all(secular(middle, p, q) < 0 for p,q in ((p_e,p_f), (p_f,p_e))):
            lo = middle
        else:
            hi = middle
    return lo, level*(p_e+p_f)


def verify(root):
    char_sum,level = character_sum_bound(root)
    assert char_sum > arb(".00058130") and level>char_sum/2
    print("任意三平面的奇夸克下界", (char_sum/2).lower().str(13))
    reflection = run_path(str(root / "ch39-fixed-reflection-quark-spectral-neighbourhood.py"))
    reference, _ = reflection["reflection_spectrum"](root)
    ranktwo = run_path(str(root / "ch39-rank-two-generator-up-quark.py"))
    centers = list(map(arb, ranktwo["CENTERS"]))
    windows = [arb(a).union(arb(b)) for a,b in ranktwo["WINDOWS"]]
    for name, masses, d_floor, s_floor in (
        ("上型中心", centers, ".00019020", ".00064244"),
        ("上型三窗", windows, ".00017876", ".00061116"),
    ):
        down, strange = bounds(reference, masses, char_sum, level)
        assert down > arb(d_floor) and strange > arb(s_floor)
        print(name, "下、奇夸克严格下界", down.lower().str(13), strange.lower().str(13))
    print("下型两字符和严格大于", char_sum.lower().str(13))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)

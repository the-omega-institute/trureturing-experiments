"""十二壳弱差全过门要求极点的去谱底权至少减少 18.69%。"""
from itertools import combinations, product
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def adj(matrix):
    return matrix.conjugate().transpose()


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    graph = lib("three-ray-chiral-boundary")
    frames, rotation = lib("cyclic-two-characters-nine")["characters"](graph["geo"]["sites"])
    lib("eleven-weak-eight-sites")["check_cyclic_geometry"](graph, rotation)
    frame = frames[0]
    native = lib("native-rh-ray-pvm")
    gram = native["gram_cholesky"](root)[3]
    _, _, finite, _, cover = lib("ninth-weak-two-line-threshold")["verify"](root)
    ctx.prec = 320
    floor, middle, top = levels
    planck, mz = arb("1.22089e19"), arb("91.1876")
    source = fifth * columns * acb_mat([[(middle-floor).sqrt(), 0], [0, (top-floor).sqrt()]])
    pole = source.tolist()[14]
    pole_norm = sum(abs(z)**2 for z in pole).sqrt()
    assert abs(frame[14,6]-1) < arb("1e-30") and pole_norm > arb(".02")
    candidates = lib("seven-cyclic-weak-horizon")["CANDIDATES"]
    centers = list(map(arb, (*lib("six-cyclic-weak-horizon")["Q"],
                        "506.579672393","504.414926074","510.694806713","520.456826617",
                        "525.342423153","535.119289427")))
    rows_by_shell, weights_by_shell = [], []
    for j in range(7,19):
        _, basis = lib("adjacent-pvm-transport")["spectrum"](graph,6,2_306_888_100_085+j)
        rows = (adj(basis)*frame).tolist()
        weights = [z.real for z in native["readout"](basis,range(19),[14,8,7,6],gram,3)]
        assert all(z.is_finite() and not z.contains(0) for z in weights)
        assert weights[0] > 0
        # 简单谱与三循环对易，每条正权谱行的非零平凡分量确定其字符。
        for row, w in zip(rows,weights,strict=True):
            if w > 0:
                assert sum(max(arb(0),abs(z).lower())**2 for z in row) > arb(".9")
        rows_by_shell.append(rows)
        weights_by_shell.append(weights)

    ref = arb("-.074256577194663365545")
    lower, near, far = arb("-1.4e-11"), arb(".020483945"), arb(".000020")
    cover((floor*planck/mz)**2, near)
    cover((arb(3000)/mz)**2, far)
    first_seven = []
    for ids, weights in zip(candidates,weights_by_shell[:7],strict=True):
        assert all(weights[i]>0 for i in ids)
        corrections = [far if i in ids else near if w>0 else lower for i,w in enumerate(weights)]
        upper = ref + sum(w*c for w,c in zip(weights,corrections,strict=True))
        assert upper < arb("-.001")
        first_seven.append(upper.upper())
    # 每壳至少一条候选低于 3000。混合选线反解极点行，和固定源行的范数矛盾。
    def radius(energy):
        return (arb(energy)/planck-floor).sqrt()

    mixed = []
    for ids in product(*candidates):
        if any(ids):
            selected = acb_mat([rows[i] for rows,i in zip(rows_by_shell[:7],ids,strict=True)])
            pole_upper = radius(3000)*sum(abs(z) for z in selected.inv().tolist()[6])
            ratio = (pole_norm/pole_upper).lower()
            assert ratio.is_finite() and ratio > 1
            mixed.append(ratio)
    assert len(mixed)==971
    print("七壳候选全高时的最大弱差、971 种混合选择的最小极点范数比",
          max(first_seven).str(12),min(mixed).str(12),flush=True)
    first = [rows[0] for rows in rows_by_shell]
    # 同一额外低线与原七壳的六条第一线配对；每个六行子集都给出必要不等式。
    other_rows = []
    for rows, weights in zip(rows_by_shell,weights_by_shell,strict=True):
        for i,w in enumerate(weights):
            if i==0 or w<0:
                continue
            bounds = []
            for six in combinations(range(7),6):
                last = acb_mat([first[h] for h in six]+[rows[i]]).inv().tolist()[6]
                bounds.append((radius(3000)*sum(abs(z) for z in last[:6]), abs(last[6])))
            assert len(bounds)==7
            other_rows.append(bounds)
    assert len(other_rows)==30

    # 六个有理复系数确定的线性组合精确消去六个自由均值，只剩固定极点。
    chosen = [1,3,4,6,8,10]
    extra = [0,2,5,7,9,11]
    seed = [("-0.3025581607248416","-0.03243356635305289"),
            ("-0.04305449548020796","0.015170175432684454"),
            ("-0.05707071653713413","-0.023777240482838333"),
            ("-0.14701253965301442","-0.11919496600313048"),
            ("-0.4313276773393039","0.000398405549657607"),
            ("0.09161180635017496","0.41935717219229424")]
    v = acb_mat([[acb(arb(a),arb(b)) for a,b in seed]])
    removed = -v*acb_mat([first[j][:6] for j in extra])*acb_mat([first[j][:6] for j in chosen]).inv()
    coeff = [acb(0)]*12
    for indices, values in ((chosen,removed.entries()),(extra,v.entries())):
        for j,c in zip(indices,values,strict=True):
            coeff[j]=c
    unavoidable = pole_norm*abs(sum(c*row[6] for c,row in zip(coeff,first,strict=True)))

    # 在第一线超过对应 cap 的同一壳，逐段包围其余所有截止到无穷高能。
    def envelope(q):
        bound = (finite((q/mz)**2)-ref).upper()*arb("1.002")+arb("1e-10")
        cover((q/mz)**2,bound)
        return bound

    def check_caps(caps, other_cap, ratio_floor, weak_ceiling):
        other_radius = radius(other_cap)
        other_ratios = []
        for seven_bounds in other_rows:
            ratios = [(pole_norm/(base+other_radius*factor)).lower() for base,factor in seven_bounds]
            assert all(v.is_finite() for v in ratios)
            other_ratios.append(max(ratios))
        assert min(other_ratios) > arb("2.98")
        available = sum(abs(c)*radius(cap) for c,cap in zip(coeff,caps,strict=True))
        dual_ratio = unavoidable/available
        norm_ratio = min(min(mixed), min(other_ratios), dual_ratio.lower())
        assert dual_ratio.is_finite() and norm_ratio > arb(ratio_floor)
        extra_high = envelope(arb(other_cap))
        uppers = []
        for q, weights in zip(caps,weights_by_shell,strict=True):
            main_high = envelope(q)
            upper = (ref+weights[0]*main_high +
                     sum(w*(extra_high if w>0 else lower) for w in weights[1:]))
            assert upper.is_finite() and upper < arb(weak_ceiling)
            uppers.append(upper.upper())
        print("额外正权线门槛 GeV、三十行排除比、十二行排除比、该壳最大弱差",
              other_cap,min(other_ratios).str(12),dual_ratio.str(13),max(uppers).str(13),flush=True)
        return dual_ratio, max(uppers), norm_ratio

    original = check_caps([q+102 for q in centers],1950,"1.0005","-.0017")
    caps = list(map(arb,(527,539,549,549,538,521,511,509,515,525,530,540)))
    _,_,required_ratio = check_caps(caps,50000,"1.10903","-.001005")
    # 若十二壳均 D >= -0.001，则极点行范数必须小于原值除以三组排除比的最小值。
    loss = 1-1/required_ratio**2
    assert loss.is_finite() and loss > arb(".1869")
    print("十二壳弱差全过门时，极点去谱底权的相对减量严格大于",loss.lower().str(12),flush=True)
    return original[:2]


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    source = run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    experiment(root,*source["verify"](root,return_operator="spectral"))

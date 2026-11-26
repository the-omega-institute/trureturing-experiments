"""用绝对归一的三个电子–μ 子手征率在减法点读出站点耦合。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import arb_mat

lepton = run_path(str(Path(__file__).with_name("ch39-polarized-lepton-rates.py")))
light, poles, q = (lepton[k] for k in ("light", "poles", "q"))
window = light["scheme"]["uncertain"]
anchor = light["scheme"]["anchor"]
physical = arb_mat([list(lepton["leptons"][a]) for a in ("A", "B")])
channels = ((0, 0), (0, 1), (1, 1))


def extract(length, z, ll, lr, rr):
    e = anchor/length
    c = poles["colour"](-e, poles["unit"])/24
    r = poles["right"](-e, poles["unit"])/4
    # 用两端正交投影后的闭式，避免区间运算重复扣除同一项。
    u_minus_v = 2*ll-3*lr+rr
    u_plus_v = 2*ll+lr-rr
    return ((27*c+z)*u_minus_v/r-u_plus_v)/(e*z)


def pair_witnesses(logs, time, reference, targets):
    (length, z, rho, _, raw_b, _), _, _ = poles["retune"]["predictions"](logs, time)
    e = anchor/length
    c = poles["colour"](-e, poles["unit"])/24
    r = poles["right"](-e, poles["unit"])/4
    side = z+9*c+r
    sr = arb_mat([[physical[i, j] for j in range(2)] for i in range(2)])
    g = arb_mat([[3*c], [-r]])
    base = sr*(arb_mat([[c, 0], [0, r]])-g*g.transpose()*(1/side))*sr.transpose()*e
    v = physical*arb_mat([[-6*c/side], [2*r/side], [1]])
    y = physical*light["mat"]["y"]
    shape, hyper = v*v.transpose(), y*y.transpose()
    left_predictions = []
    for target in targets:
        k = 1/(poles["offset"]+target/(15*length))
        for pair in combinations(channels, 2):
            design = arb_mat([[shape[i, j], hyper[i, j]] for i, j in pair])
            assert abs(design.det()) > 0
            rhs = arb_mat([[reference[i, j]-k*base[i, j]] for i, j in pair])
            strength, shift = (design.inv()*rhs).entries()
            assert strength > 0 and 0 < shift+k/24 < 1
            b = raw_b/rho
            d = (k*e/strength-4/side-b)/e
            assert d > 0
            full = light["mat"]["massive_exchange"](e, z, b, d)*k + light["mat"]["hypercharge"]*(shift+k/24)
            assert light["mat"]["positive"](full)
            trial = physical*full*physical.transpose()
            assert all(trial[i, j].overlaps(reference[i, j]) for i, j in pair)
            if pair == channels[1:]:
                left_predictions.append(trial[0, 0])
    return left_predictions


def verify():
    assert poles["verify_ratio_monotonicity"]() > 0
    time = poles["study"]["nominal"][2]
    targets = [poles["study"]["comp"]["model"]["target"](n, time) for n in (11430, 11431, 11432)]
    for case, logs in enumerate(poles["retune"]["log_sets"]):
        graph = light["graphs"](logs, time)[0]
        m1, m2 = poles["reference_readings"](logs)[:2]
        absolute = physical*graph*physical.transpose()
        pair_witnesses(logs, time, absolute, (targets[0], targets[2]))
        diff = window(q("1e-10"))
        length, z, _ = poles["charged_geometry"](m1*diff, m2*diff)
        magnitude = [(absolute[i, j]**2*window(q("1e-9"))).sqrt() for i, j in channels]
        a, b, c = magnitude
        for sign in (1, -1):
            k = extract(length, z, a, sign*b, c)*window(q("2e-5")).sqrt()
            target = 15*length*(1/k-poles["offset"])
            assert k > 0
            if sign == 1:
                assert targets[2] < target < targets[0] and target.contains(targets[1])
                assert q("3.2939705") < target < q("3.2941822")
            else:
                assert target < 0
            print("case", case, "LR sign", sign, "target endpoints", target.lower().str(12),
                  target.upper().str(12), flush=True)


if __name__ == "__main__":
    verify()

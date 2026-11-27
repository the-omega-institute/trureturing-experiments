"""相对效率极大极小权：在弦度量锥 c_s=(1,1,s)，s in [0,1] 上，e(w)=min_s g_s(w)/G(s) 的最大值 e* 的两侧界。

沿用 ch39-margin-metric-family-vacuum.py 的对象、覆盖与行（本脚本从同名 .json 的 "maximin" 段读证书）。
  1. 下界：有理见证 W_e；把 [0,1] 分成有理单元 [a,b]，单元上 G(s) <= G(a) <= U_a、sG(s) <= bG(b) <= bU_b，
     故 G(s) min_C rho_{C,s}^2 <= min_C(U_a(r1_C+r2_C) + bU_b r3_C)；全环面覆盖证明
     N_{W_e} >= min_C(gam(r1_C+r2_C) + sig r3_C)，gam >= e_lo U_a，sig >= e_lo b U_b，于是 e(W_e) >= e_lo。
     每个 U_a 由对偶证书给出。
  2. 上界：对若干节点 s_k，节点见证的覆盖给出 G(s_k) >= L_k；e(w) >= e 推出 a_j w >= e L_k min_C(R12_C + s_k R3_C)，
     对偶证书给出 e* <= e_hi < 1。
  3. 一切 e(w) >= e_lo 的权：w_105、w_104 的区间；由 g_0(w) >= e_lo L_0 > 1/1000 与漂移证书，它们都是严格共同真空。
"""
from pathlib import Path
from runpy import run_path
import json
import sys

from flint import arb

HERE = "ch39-margin-metric-family-vacuum"


def verify(root):
    fam = run_path(str(root/(HERE+".py")))
    base = run_path(str(root/"ch39-quadratic-growth-margin-vacuum.py"))
    model = base["Model"](root)
    hessian = base["hessian"]
    alldata = json.loads((root/(HERE+".json")).read_text())
    data = alldata["maximin"]
    model.set_phase(arb(alldata["t"]))
    rows = fam["Rows"](model, data["rows"])
    frac = fam["frac"]
    J = range(len(rows.data))

    # 上界 U(s)：对一切 w in P，g_s(w) <= U(s)。
    U = {}
    for key, cert in data["uppers"].items():
        s = frac(key)
        numer, norm = model.dual_value(cert, rows.level([(j, arb(1), s) for j in J]), {}, arb(0))
        assert norm > 0
        U[key] = numer/norm
    # 见证 W_e 与单元覆盖。
    e_lo = arb(data["e_lo"])
    W = model.complete({int(q): arb(v) for q, v in data["free"].items()}, data["basis"])
    modes = [(q, [W[q]*c for c in model.bare[q]]) for q in W]
    cells = data["cells"]
    assert cells[0]["a"] == "0/1" and cells[-1]["b"] == "1/1"
    assert all(x["b"] == y["a"] for x, y in zip(cells, cells[1:]))
    for c in cells:
        a, b = frac(c["a"]), frac(c["b"])
        assert a < b
        gam, sig = arb(c["gam"]), arb(c["sig"])
        assert gam >= e_lo*U[c["a"]] and sig >= e_lo*b*U[c["b"]]
        fam["cover"](hessian, modes, model.x, model.y, gam, sig)
    assert data["e_lo"] == "0.916" and len(cells) == 40 and len(W) == 13
    assert abs(W[105]-arb("0.11422969")) < arb("1e-8") and abs(W[104]-arb("0.02558634")) < arb("1e-8")
    assert W[105] > arb("0.112967")     # 高于全部 21 个节点最优集的 w_105
    print("相对效率下界：e(W_e) >= e_lo、单元数、W_e 支撑、w_105、w_104", e_lo.str(8), len(cells), len(W),
          W[105].str(8), W[104].str(8), flush=True)
    # 节点下界 L_k（覆盖）与 e* 上界。
    L = {}
    for nd in data["nodes"]:
        s = frac(nd["s"])
        Wk = model.complete({int(q): arb(v) for q, v in nd["free"].items()}, nd["basis"])
        mk = [(q, [Wk[q]*c for c in model.bare[q]]) for q in Wk]
        Lk = arb(nd["L"])
        fam["cover"](hessian, mk, model.x, model.y, Lk, Lk*s)
        L[nd["s"]] = (s, Lk)
    pairs = [(j, L[k][1], L[k][1]*L[k][0]) for k, j in data["pairs"]]
    lv = rows.level(pairs)
    numer, norm = model.dual_value(data["upper"], lv, {}, arb(0))
    assert norm > 0
    e_hi = numer/norm
    assert e_lo < e_hi < arb("0.940879")
    assert sorted(L) == sorted(f"{k}/20" for k in (0, 2, 3, 9, 10, 11, 12, 14, 20))
    print("相对效率上界 e* <= e_hi", e_hi.str(10), flush=True)
    # e(w) >= e_lo 的权：区间。
    lv = rows.level([(j, e_lo*al, e_lo*be) for j, al, be in pairs])
    enc = {}
    for key, phi in (("w105", {105: arb(1)}), ("-w105", {105: arb(-1)}), ("w104", {104: arb(1)}), ("-w104", {104: arb(-1)}),
                     ("drift", model.mu)):
        enc[key] = model.dual_value(data["enclose"][key], lv, phi, arb(1))[0]
    assert -enc["-w105"] > arb("0.110406") and enc["w105"] < arb("0.120795")
    assert -enc["-w104"] > arb("0.0242976") and enc["w104"] < arb("0.0279012")
    assert enc["drift"] < arb("-2.0763")
    # 严格性：g_0(w) >= e_lo L_0 > 1/1000。
    assert data["nodes"][0]["s"] == "0/20" and data["nodes"][0]["L"] == "6.348723"
    assert e_lo*L["0/20"][1] > arb(1)/1000
    print("e(w) >= e_lo 的权：w_105、w_104 区间与漂移上界", (-enc["-w105"]).str(7), enc["w105"].str(7),
          (-enc["-w104"]).str(7), enc["w104"].str(7), enc["drift"].str(6))
    return dict(e_lo=e_lo, e_hi=e_hi, enc=enc, W=W)


if __name__ == "__main__":
    here = Path(__file__).resolve().parent
    verify(Path(sys.argv[1]) if len(sys.argv) > 1 else here)

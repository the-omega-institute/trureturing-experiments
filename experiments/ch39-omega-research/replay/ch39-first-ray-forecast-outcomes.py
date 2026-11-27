"""用实际十九阶磁拉普拉斯的三个本征向量读取预登记的五行二十五项。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb


def verify(root):
    plan = run_path(str(root/"ch39-first-ray-forecast-phases.py"))
    graph = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    spectrum = run_path(str(root/"ch39-adjacent-pvm-transport.py"))["spectrum"]
    sites = graph["boundary_sites"]()
    assert set(sites) == {6,7,8}
    hits = 0
    others_lower = arb(1)
    predicted_upper = arb(0)
    for predicted,number in plan["FORECAST"]:
        weak = []
        bounds = []
        for charge,_ in plan["FORECAST"]:
            levels, basis = spectrum(graph,charge,number)
            assert levels[2] < 9 < levels[3]
            rows = acb_mat([[basis[i,j] for j in range(3)] for i in sites])
            values = (rows.conjugate().transpose()*rows).eig(algorithm="rump")
            assert all(v.is_finite() and v.imag.contains(0) and v.real > 0 for v in values)
            lo = min(v.real.sqrt().lower() for v in values)
            hi = min(v.real.sqrt().upper() for v in values)
            if hi < arb(".001"):
                weak.append(charge)
            else:
                assert lo > arb(".001")
            bounds.append((lo,hi))
            if charge == predicted:
                predicted_upper = max(predicted_upper,hi)
            else:
                others_lower = min(others_lower,lo)
        hits += weak == [predicted]
        print("编号、预测荷、实测读弱荷、按五荷顺序的奇异值上端",
              number, predicted, weak, *(hi.str(12) for _,hi in bounds), flush=True)
    print("成功定位唯一读弱荷的行数", hits, "共", len(plan["FORECAST"]))
    assert hits == 5 and predicted_upper < arb(".000315") and others_lower > arb(".03384")
    return hits


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)

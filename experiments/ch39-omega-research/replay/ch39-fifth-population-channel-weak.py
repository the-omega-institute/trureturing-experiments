"""保持第五壳全部原生谱占据数的量子信道，在三个后继测试壳逐壳失配弱差。"""
from pathlib import Path
from runpy import run_path
from flint import arb


def adj(m):
    return m.conjugate().transpose()


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    weak = lib("three-shell-rank-one-obstruction")
    # 此上界计算只用源的最低、最高振幅和原生实权。
    # 对每个壳，若其所有候选均不低于所列能标，完整弱核的和严格为负。
    weak["experiment"](root, fifth, columns, levels)
    calculate = lib("adjacent-pvm-transport")["spectrum"]
    graph = lib("three-ray-chiral-boundary")
    floor, heavy, top = levels
    planck = arb("1.22089e19")
    for j, ids in enumerate(weak["CANDIDATES"]):
        _, target = calculate(graph, 6, 2_306_888_100_092+j)
        overlap = adj(target)*fifth
        bounds = []
        for i in ids:
            terms = [[abs(overlap[i,k]*columns[k,h]) for k in range(19)] for h in (0,1)]
            # Kraus 环境的第 k 个向量范数为 1。每个重源列独立应用反三角不等式。
            gain = []
            for row in terms:
                gain.append(max([arb(0)]+[(v-sum(w for k,w in enumerate(row) if k!=l)).lower()
                                          for l,v in enumerate(row)]))
            bounds.append((floor + (heavy-floor)*gain[0]**2 + (top-floor)*gain[1]**2).lower())
        energy = planck*min(bounds)
        print("壳号、原生候选各线的统一能量下界（GeV）", j+7, energy.lower().str(15), flush=True)
        assert energy > arb("1e15") > weak["THRESHOLDS"][j]


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    source = run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    experiment(root, *source["verify"](root, return_operator="spectral"))

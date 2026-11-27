"""依次保持第五、第六壳原生占据数的两次信道不能给出第七至九壳任一弱端点。"""
from pathlib import Path
from runpy import run_path
from flint import arb


def adj(m):
    return m.conjugate().transpose()


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    weak = lib("three-shell-rank-one-obstruction")
    weak["experiment"](root, fifth, columns, levels)
    native = lib("adjacent-pvm-transport")["spectrum"]
    graph = lib("three-ray-chiral-boundary")
    sixth = native(graph, 6, 2_306_888_100_091)[1]
    between = adj(sixth)*fifth
    floor, heavy, top = levels
    for j, candidates in enumerate(weak["CANDIDATES"]):
        target = native(graph, 6, 2_306_888_100_092+j)[1]
        left = adj(target)*sixth
        bounds = []
        for line in candidates:
            amplitudes = []
            for h in range(2):
                paths = [abs(left[line,b]*between[b,a]*columns[a,h])
                         for a in range(19) for b in range(19)]
                largest = max(range(len(paths)), key=lambda k: paths[k].lower())
                lower = paths[largest]-sum(z for k,z in enumerate(paths) if k!=largest)
                amplitudes.append(max(arb(0),lower.lower()))
            bounds.append((floor + (heavy-floor)*amplitudes[0]**2
                           + (top-floor)*amplitudes[1]**2).lower())
        energy = arb("1.22089e19")*min(bounds)
        print("先第五壳后第六壳的保占据信道、目标壳、所有候选截止下界（GeV）",
              j+7, energy.lower().str(15), flush=True)
        assert energy > arb("1e14") > weak["THRESHOLDS"][j]


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    source = run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    experiment(root, *source["verify"](root, return_operator="spectral"))

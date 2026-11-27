"""右外壳上型完整匹配要求七维三循环不变空间向其十二维正交补转移。"""
from pathlib import Path
from runpy import run_path

from flint import arb


def native_bounds(root):
    prior = run_path(str(root / "ch39-character-preserving-up-ckm-obstruction.py"))
    prior["transport_bounds"](root)
    trained = run_path(str(root / "ch39-common-embedding-central-right-ckm.py"))
    native = run_path(str(root / "ch39-free-four-site-native-positive-spectrum.py"))["native"]
    n = 2_306_888_100_086
    _, labels, _ = trained["verify_nodes"](root, ((n, ".0294"),), ".1", (),
                                          ((6, 7, 10, 13), (6, 7, 10, 12)))
    rows = native(root, n)[0]
    triv = [j for j,r in enumerate(rows) if (r[0,0]+r[1,1]).real > 0]
    assert len(triv) == 7
    w = labels[n,-4]
    assert all(0 < w[j] < arb(".412") for j in triv)
    assert all(arb(".00191969") < w[j] < arb(".947555") for j in range(19) if j not in triv)


def probability_ratios(m):
    u, c, t = m
    e0, f1 = arb(".00191969"), arb(".947555")
    h0, h1 = arb(".0077")-u, arb(".412")
    # 从多项式在三个质量处的取值比界住两个概率比，分母全有固定正下界。
    r2 = (c-e0)/(e0-u) * u/(c-u) * (h0-u)/(h0-c) * (t-c)/(t-u)
    r3 = (t-e0)/(e0-u) * (f1-t)/(f1-u) * u*h1/((t-u)*(t-h1)) * (t-c)/(c-u)
    assert r2 < arb(".00206") and r3 < arb(".07374")
    print("整个字符盒给出的两个弱轴概率比上界", r2.str(15), r3.str(15), flush=True)
    return {1: arb(".00206"), 2: arb(".07374")}


def verify(root):
    native_bounds(root)
    model = run_path(str(root / "ch39-triple-geometry-pencil.py"))
    lib, geom, mixing, mass, up, down, _ = model["setup"](root)
    ratios = probability_ratios(mass[0])
    normalized = [v/mass[0][-1] for v in mass[0]]
    for chart, (a,b) in enumerate(((up,down), (down,up))):
        for branch in range(3):
            stack = [(arb(-1), arb(1), 0)]
            leaves = cancelled = deepest = 0
            while stack:
                lo, hi, depth = stack.pop()
                lines = lib["lines"](a, b, lo, hi, (branch,), geom)
                if branch in lines:
                    y = lines[branch]
                    p = [abs(v)**2 for v in (mixing*y).entries()]
                    done = any(p[0] > ratio*p[j] for j,ratio in ratios.items())
                    if not done and chart == 0:
                        # 特别是 t=0 附近，消去两个非主分量所共有的 t，以其极限包围。
                        dy = [abs(v)**2 for v in (mixing*(down*y)).entries()]
                        radius = max(abs(lo),abs(hi)).upper()
                        eigen = normalized[branch] + arb(0,radius)
                        delta = [(v-eigen)**2 for v in normalized]
                        done = any(dy[0]*delta[j] > ratio*dy[j]*delta[0]
                                   for j,ratio in ratios.items())
                        if done:
                            cancelled += 1
                    if done:
                        leaves += 1
                        deepest = max(deepest, depth)
                        continue
                assert depth < 26, (chart, branch, lo, hi)
                mid = (lo+hi)/2
                stack.extend(((lo,mid,depth+1), (mid,hi,depth+1)))
            print("铅笔图、升序支、叶区间、消去参数的叶区间、最大深度",
                  chart, branch, leaves, cancelled, deepest, flush=True)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)

"""两部门完整边相位中，保留中心星时的长边最大角误差。"""
from itertools import permutations, product
from pathlib import Path
from runpy import run_path

from flint import acb, arb


PHASE_TRILLIONTHS = (-589008033029,-1374438668050,-634060530989,57991527584,
             -1236007593642,57991527584,-3648161005871,-805905093363,
             -158905532750,57991527584,57991527584,57991527584,57991527584,
             -351904688965,0,57991527584,57991527584,-521011918124,57991527584)


def verify():
    shell = run_path(str(Path(__file__).with_name("ch39-arc-relative-determinant.py")))
    heights = [acb.zeta_zero(n).imag for n in (2345666,964890)]
    targets, clocks = {}, {}
    for i,j,d,z in shell["edges"]:
        if d<=372:
            continue
        phase = acb(1) if 14 in (i,j) else z**-6
        if d in (724,804):
            omega = shell["theta"](d)-shell["theta"](244)
            phase *= acb(0,omega*(-4*heights[0]-2*heights[1])).exp()
            clocks[i,j] = phase.arg()
        targets[i,j] = phase.arg()
    assert all(targets[i,14]==0 for i in range(6))
    pi = arb.pi()
    cap = arb("1.39")
    assert 4*cap < 2*pi
    lifts = []
    for j in (6,7,8):
        edges = [edge for edge in clocks if edge[1]==j]
        assert len(edges)==4
        patterns = []
        for turns in product((-1,0,1),repeat=3):
            lifted = {e:clocks[e]+2*pi*k for e,k in zip(edges,(0,*turns),strict=True)}
            if any(x-y>4*cap for x,y in permutations(lifted.values(),2)):
                continue
            patterns.append(lifted)
        lifts.append(patterns)

    nodes = (6,7,8,14)
    cycles = [seq for length in (2,3,4) for seq in permutations(nodes,length)
              if seq[0]==min(seq)]
    threshold = arb("1.37443866803")
    for triple in product(*lifts):
        clock_lift = {e:t for pattern in triple for e,t in pattern.items()}
        # 对每个外点消去同一实角，得 x_k−x_j ≤ offset[j,k]+2r。
        offset = {}
        for i in range(6):
            labels = {j:t for (site,j),t in clock_lift.items() if site==i}
            labels[14] = arb(0)
            for j,k in permutations(labels,2):
                cost = labels[k]-labels[j]
                offset[j,k] = offset.get((j,k),cost).min(cost)
        assert any(-sum(offset[j,k] for j,k in zip(seq,seq[1:]+seq[:1],strict=True))
                   > 2*len(seq)*threshold for seq in cycles)

    angles = [arb(v)/10**12 for v in PHASE_TRILLIONTHS]
    upper = arb("1.37443866806")
    for (i,j),t in targets.items():
        assert abs(acb(0,t+angles[i]-angles[j]).exp().arg()) < upper
    print("完整背景、中心星保留时的最优长边角误差",threshold,upper)


if __name__ == "__main__":
    verify()

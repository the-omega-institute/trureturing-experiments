"""独立轻线只能调相时，达到参考 CKM 混合至少需要几个单线酉相位门。"""
from itertools import permutations
from pathlib import Path
from runpy import run_path

from flint import arb


def verify():
    root = Path(__file__).resolve().parent
    model = run_path(str(root/"ch39-shell-clock-ckm-blocks.py"))
    for title, target in (("参考中心",model["standard_matrix"](*model["centers"])),
                          ("整个参考盒",model["target_box"]())):
        distances = []
        for perm in permutations(range(3)):
            # 同时重排两轴保住对称性，故只需穷举一侧的六个相对匹配。
            skew = [(abs(target[i][perm[j]])-abs(target[j][perm[i]]))
                    for i in range(3) for j in range(i)]
            bound = (sum(v*v for v in skew)/2).sqrt()
            distances.append(bound)
        best = min(range(6),key=lambda i:float(distances[i]))
        assert all(distances[best]<d for i,d in enumerate(distances) if i!=best)
        print(title,"单一秩一酉相位门的误差下界及最佳相对匹配",
              distances[best].str(17),list(permutations(range(3)))[best])
        assert distances[best] > arb(".0024")


if __name__=="__main__":
    verify()

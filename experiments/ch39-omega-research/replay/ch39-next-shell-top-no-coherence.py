"""完全去相干时下一壳的最高上型质量，以及对任意相干比例通用的固定输运。"""
from itertools import permutations
from pathlib import Path
from runpy import run_path

from flint import acb, arb


def verify(root):
    memory = run_path(str(root/"ch39-next-shell-partial-coherence-up-threshold.py"))
    _,diagonal,target,error = memory["source_data"](root)
    entries,_ = memory["target_data"](target)
    small_eigen = run_path(str(root/"ch39-next-shell-arbitrary-up-transport-floor.py"))["small_eigen"]
    weights = []
    for j in range(19):
        center = diagonal[j,j].real.sqrt()
        assert center>error
        weights.append(((center-error)**2).union((center+error)**2))
    lower = sorted((v.lower() for v in weights),reverse=True)
    upper = sorted((v.upper() for v in weights),reverse=True)
    # 最大本征值等于负矩阵最小本征值的相反数；只有七条目的谱线贡献平凡块。
    def max_trivial(w):
        return max((-small_eigen(tuple(-v for v in perm),entries)).upper() for perm in permutations(w[:7]))
    trivial_upper = max_trivial(upper)
    root3 = arb(3).sqrt()
    omega = acb(-1,root3)/2
    strengths,limits = [],[]
    for o in (omega,omega.conjugate()):
        s = [abs(target[8,j]+o*target[7,j]+o**2*target[6,j])**2/3 for j in range(19)]
        assert sum(v>0 for v in s)==6 and sum(s).contains(1)
        lo = sorted((max(arb(0),v.lower()) for v in s),reverse=True)
        hi = sorted((v.upper() for v in s),reverse=True)
        strengths.append(s)
        limits.append(sum(a*b for a,b in zip(lower,lo,strict=True)).union(
                      sum(a*b for a,b in zip(upper,hi,strict=True))))
    hmax,second = limits
    assert trivial_upper<arb(".909321") and second<arb(".917832")
    assert max(trivial_upper,second.upper())<hmax
    assert arb(".92832724")<hmax<arb(".92832735")
    print("完全去相干下的上型上端，全局包围和另外两块上界",
          hmax.lower().str(14),hmax.upper().str(14),second.upper().str(14),trivial_upper.str(14),flush=True)

    # 源、目的层均按一基编号列出。每对中最大的概率配最大的源振幅。
    src = [j-1 for j in (17,14,2,4,5,8)]
    dest = [j-1 for j in (17,14,2,5,8,11)]
    others_src = [i for i in range(19) if i not in src]
    others_dest = [i for i in range(19) if i not in dest]
    assert all(weights[i]>weights[j] for i,j in zip(src,src[1:]))
    leftover = max(weights[j].upper() for j in others_src)
    assert weights[src[-1]]>leftover and leftover<arb(".408")
    s0,s1 = strengths
    assert all(s0[i]>s0[j] for i,j in zip(dest,dest[1:]))
    assert all(s0[j].contains(0) for j in others_dest)
    assert all(s1[j].contains(0) for j in dest)
    assign = dict(zip(dest+others_dest,src+others_src,strict=True))
    assert len(set(assign.values()))==len(assign)==19
    y = [weights[assign[j]] for j in range(19)]
    high,low = [sum(a*b for a,b in zip(row,y,strict=True)) for row in (s0,s1)]
    assert high.overlaps(hmax) and high>arb(".8659") and low<leftover
    goal = arb(".8616")
    p = (goal-low)/(high-low)
    assert 0<p<1
    numerator = 8785346
    fixed = arb(numerator)/10**7
    fixed_top = fixed*high+(1-fixed)*low
    assert abs(fixed_top-goal)<arb("1e-7")
    print("固定谱基排列的另一字符、余十三层上界、固定混合概率分子和顶质量",
          low.str(12),leftover.str(12),numerator,fixed_top.lower().str(14),
          fixed_top.upper().str(14),flush=True)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)

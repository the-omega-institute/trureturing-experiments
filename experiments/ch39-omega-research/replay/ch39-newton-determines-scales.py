"""以牛顿矩逐背景消去独立标尺，允许任意逐线连续位置概率分布。"""
from pathlib import Path
from runpy import run_path

from flint import arb

ROOT = Path(__file__).resolve().parent


def check_model():
    shared = run_path(str(ROOT / 'ch39-shared-twelve-node-gauge-gap.py'))
    shell = run_path(str(ROOT / 'ch39-native-rh-ray-pvm.py'))['shell_weights']
    rebuild = run_path(str(ROOT / 'ch39-two-zero-cutoff-gravity.py'))['gravity_fields']
    _, lib, backgrounds = rebuild(ROOT, shared['HEAVY'], tight_neutral=True)
    constants = []
    shared_positions = None
    for number, rows in zip(shared['SHELL'], backgrounds, strict=True):
        edges = shared['paired_lines'](ROOT, shell, number)
        assert all(0 < left and (left < right or left == right) and right <= 1
                   for left, right, _ in edges)
        # E₁(z/x) 随 x 递增。去掉零质量行，其齐次误差恒等于零。
        bound = sum(abs(weight) * sum(abs(alpha) * z * lib['expint'](1, z / right)
                    for _, alpha, z in rows[1:]) for _, right, weight in edges) / (2*arb.pi())
        assert rows[0][2].is_zero()
        constants.append(bound.upper())
        if shared_positions is None:
            assert all(not w.contains(0) for _, _, w in edges)
            shared_positions = [w < 0 for _, _, w in edges]
        # 第一背景负权线取均值端点、正权线取乘积端点；两个背景共用选择表。
        locations = [(right if mode else left, weight)
                     for (left, right, weight), mode in zip(edges, shared_positions, strict=True)]
        def read(c):
            return sum(weight*c*node*alpha*lib['expint'](2,z/(c*node))
                       for node,weight in locations for _,alpha,z in rows)/(2*arb.pi())
        assert read(arb(1)) > 500
        lo,hi=arb('0.001'),arb('0.003')
        assert read(lo)<1<read(hi)
        # 插值求种子，最后只用经区间重新计算的两端点来夹根。
        for _ in range(4):
            a,b=read(lo),read(hi)
            guess=(lo+(1-a)*(hi-lo)/(b-a)).mid()
            if read(guess)<1:
                lo=guess
            else:
                hi=guess
        guard=arb('1e-9')
        lo,hi=guess-guard,guess+guard
        assert read(lo)<1<read(hi)
        print('同一端点选择的背景、唯一牛顿标尺所在区间',number,lo.union(hi).str(13),flush=True)
    assert constants[0] < arb('0.241') and constants[1] < arb('0.289')
    print('两背景全部标尺、全部逐线分布的 |cK′−K| 上界',
          *[value.str(13) for value in constants])
    print('在 K=1 时对 log c 的导数下界',
          *[(1-value).str(13) for value in constants])


if __name__=='__main__':
    check_model()

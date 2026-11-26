"""用有理源和输入核对弱三分量局部动能的独立方向。"""
from fractions import Fraction as Q


def verify():
    # 中性 Higgs 的期望荷依次对应 Q,T₃,T₁,T₂。
    expectation = [Q(0), Q(-1, 2), Q(0), Q(0)]
    assert [[a*b for b in expectation] for a in expectation] == [
        [0,0,0,0], [0,Q(1,4),0,0], [0,0,0,0], [0,0,0,0]]
    v, mz, planck = map(Q, ("246.22", "91.1876", "1.22089e19"))
    rho = (v/planck)**2/2
    def shift(energy):
        return rho**2*(energy**2-mz**2)/(4*planck**2)
    assert shift(mz) == 0
    for energy, lower, upper in (
        ("1e13","6.936e-81","6.937e-81"),
        ("1.45137757e14","1.461e-78","1.462e-78"),
        ("1e17","6.936e-73","6.937e-73"),
    ):
        value = shift(Q(energy))
        assert Q(lower) < value < Q(upper)
        print(energy, "GeV，单位系数给 A₃₃ 的位移", float(value))


if __name__ == "__main__":
    verify()

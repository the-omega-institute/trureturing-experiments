"""右外壳的全部九站支撑均排除轻上夸克参考窗。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    exceptions = [tuple(sorted((0,1,2,3,4,5,14) + pair)) for pair in combinations((11,12,18),2)]
    assert len(exceptions) == 3
    check = run_path(str(root/"ch39-up-eight-site-cover.py"))["verify"]
    check(root, 9, "ch39-up-nine-site-cover.tsv", ".00000629", exceptions)
    special = run_path(str(root/"ch39-up-nine-site-coherent-pair.py"))
    assert special["SUPPORT"] in exceptions
    # 三个剩余支撑在三循环下彼此等价，两端算子、源振幅和四站子空间协变。
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    sites = shell["geo"]["sites"]
    def rotate(support):
        return tuple(sorted(sites.index((sites[i][1],sites[i][2],sites[i][0])) for i in support))
    one = special["SUPPORT"]
    assert {one,rotate(one),rotate(rotate(one))} == set(exceptions)
    assert set(rotate((14,8,7,6))) == {14,8,7,6}
    special["verify"](root)
    print("全部 92378 个九站支撑已排除，块外恒等时至少需要十站。")

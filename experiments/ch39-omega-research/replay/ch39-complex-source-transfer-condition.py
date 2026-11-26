"""从四种实际算术三源到单位有界的复味空间时，求等变可逆输运的最小条件数。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify(root):
    get_gram = run_path(str(root / "ch39-native-rh-ray-pvm.py"))["gram_cholesky"]
    normal = run_path(str(root / "ch39-four-gram-triplet-ckm.py"))["normal"]
    ctx.prec = 320
    table = get_gram(root)
    a = (arb(".00384")*arb(".0384")).sqrt()
    p = arb(10)/11
    b = (1-a*a).sqrt()
    target = [a,b*p.sqrt(),b*(1-p).sqrt()]
    for omit in range(4):
        n = [v.real for v in normal(table,omit)]
        fixed = (n[0]**2+sum(n[1:])**2/3).sqrt()
        rotating = ((1-fixed**2)/2).sqrt()
        gains = [u/v for u,v in zip((fixed,rotating,rotating),target,strict=True)]
        assert gains[1] < gains[2] < gains[0]
        condition = gains[0]/gains[1]
        assert condition > 60
        assert (gains[2]/gains[1]).overlaps(arb(10).sqrt())
        print("省略源、最小二范数条件数、两个字符的增益比",omit,condition.str(14),
              (gains[2]/gains[1]).str(10))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)

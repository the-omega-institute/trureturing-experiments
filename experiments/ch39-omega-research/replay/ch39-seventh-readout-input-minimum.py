"""在同一组实际零点上比较两份产生相同六个43位输入码的正权表。"""
from pathlib import Path
from runpy import run_path
from flint import arb, ctx


def verify():
    root = Path(__file__).parent
    dictionary = run_path(str(root / "ch39-seventh-time-zero-dictionary.py"))
    radial = run_path(str(root / "ch39-six-time-neutral-normalization.py"))
    actual = run_path(str(root / "ch39-seven-actual-neutral-times.py"))
    support = dictionary["upper_support"]
    coordinates = {n: x/(1+x) for n, (_, x, _) in
                   zip(support, actual["shell_nodes"](2_903_600_000, support), strict=True)}
    ctx.prec = 576
    times = list(radial["neutral"]["model"]["shared"]["clock"]["clock"]["times"]()[:7])
    times[6] += arb("1.360900010e-8")
    scale = 2**43
    codes = [(t*scale).floor().unique_fmpz() for t in times[1:]]
    assert all(c is not None for c in codes)
    inset = arb(1)/2**20
    perturbed = [arb(1), *[(n+(1-inset if j%2==0 else inset))/scale
                            for j,n in enumerate(codes)]]
    for values, output in ((times, 32071996606), (perturbed, 32071996605)):
        assert [(t*scale).floor().unique_fmpz() for t in values[1:]] == codes
        den, moments = radial["compact_moments"](values)
        reading, _ = dictionary["certified_rule"](support, coordinates, values, moments, den)
        assert isinstance(reading, arb)
        assert (reading*2**35).floor().unique_fmpz() == output
        print("同七枚正权零点的43位输入与35位输出", *codes, output)


if __name__ == "__main__":
    verify()

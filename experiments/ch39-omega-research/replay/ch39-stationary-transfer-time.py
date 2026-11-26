"""Concentration and first-passage obstructions for finite stationary power-weight models."""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

uphill = run_path(str(Path(__file__).with_name("ch39-uphill-transfer-radius.py")))


def verify():
    ctx.prec = 144
    scores, _ = uphill["scores"]()
    target = uphill["TARGET"]
    ideal = {n: x for n, (x, _) in scores.items()}
    robust = {n: x for n, (_, x) in scores.items()}
    odds_limit = arb(1)/99
    assert (ideal[10871]/ideal[target])**5967 > odds_limit
    odds = sum((w/ideal[target])**5968 for n, w in ideal.items() if n != target)
    assert odds < odds_limit
    exponent = arb(99).log() / (robust[target].upper()/robust[10871].lower()).log()
    assert exponent > arb("5616.39")
    odds = sum((w.upper()/robust[target].lower())**6367 for n, w in robust.items() if n != target)
    assert odds < odds_limit
    print("Ideal least integer exponent 5968; finite filter: necessary > 5616.39, sufficient 6367")

    # Every index path from 11499 to 10875 with these step bounds must visit the
    # specified band. Its entries lie in the guaranteed feasible interval.
    for radius, left, decimal_floor in [(1,11492,6423), (2,11351,5801), (3,11480,3803),
                                         (8,11275,209), (16,11464,30), (32,11442,12)]:
        cap = max(robust[n].upper() for n in range(left, left+radius))
        ratio = robust[11499].lower()/cap
        log_time = (exponent*ratio.log()-(2*arb(radius)).log())/arb(10).log()
        assert log_time > decimal_floor
        print("Jump limit", radius, "median time exceeds 10^", decimal_floor)

    count = uphill["LAST"]-uphill["FIRST"]+1
    assert count == 1238
    print("Radius 58 uphill-tree uniform mean first-passage bound", 116*count*(count-1)//2)


if __name__ == "__main__":
    verify()

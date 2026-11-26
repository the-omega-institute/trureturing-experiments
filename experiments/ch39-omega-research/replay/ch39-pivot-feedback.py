"""Stationary zero-label concentration under one filtered reference reflection."""
from pathlib import Path
from runpy import run_path
from itertools import islice

from flint import acb, arb

pulses = run_path(str(Path(__file__).with_name("ch39-finite-port-pulses.py")))


def verify():
    scores = {}
    weights = {}
    for index, zero in enumerate(acb.zeta_zeros(10815, 1278), 10815):
        series = pulses["pulse_records"](pulses["generator"](zero.imag))
        scores[index] = next(islice(series, 6, 7))[0]
        weights[index] = 1/(arb(1)/4 + zero.imag**2)
    order = sorted(scores, key=lambda n: float(scores[n].mid()))
    assert all(scores[i] < scores[j] for i, j in zip(order, order[1:]))
    runner, winner = order[-2:]
    assert winner == 11431 and 0 < scores[order[0]] and scores[winner] < 1
    ratios = [(1-scores[winner])/(1-scores[n]) for n in order[:-1]]

    def stationary_share(repetitions):
        return 1/(1+sum(q**repetitions for q in ratios))

    lo, hi = 0, 32768
    target = arb(99)/100
    assert stationary_share(hi) > target
    while hi-lo > 1:
        mid = (lo+hi)//2
        probability = stationary_share(mid)
        if probability > target:
            hi = mid
        else:
            assert probability < target
            lo = mid
    assert stationary_share(lo) < target < stationary_share(hi)
    assert hi == 10175
    escape = (1-scores[runner])**hi
    assert escape < arb("1.31e-775")
    print("one-pulse stationary share", stationary_share(1).str(12), "runner", runner)
    print("first batch above 99%:", hi, stationary_share(lo).str(12), stationary_share(hi).str(12))
    print("runner all-miss probability per batch", escape.str(9))
    print("necessary log10 batches for 90% from runner", ((arb(".9")/escape).log()/arb(10).log()).str(9))
    odds = {i: p/(1-p) for i, p in scores.items()}
    distance = 2048
    duel_share = 1/sum((x/odds[winner])**distance for x in odds.values())
    total_weight = sum(weights.values())
    w = weights[winner]/total_weight
    common_mass = min(2*w*min(weights.values())/total_weight, (1-2*w)**2)
    rounds = 2**23
    mixing_error = (-rounds*common_mass).exp()
    jump_floor = 2*scores[order[0]]*(1-scores[winner])
    # A score changes by one with probability at least jump_floor per shot pair.
    block = 2**30
    assert 8*distance**2 / jump_floor < block
    truncation = arb(rounds)/2**48
    assert duel_share-mixing_error-truncation > target
    print("pairwise pivot proposals:", duel_share.str(12), "equilibrium, common mass", common_mass.str(9))
    print("finite register lower occupancy", (duel_share-mixing_error-truncation).str(12),
          "jump floor",jump_floor.str(8))


if __name__ == "__main__":
    verify()

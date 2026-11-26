"""Bound the global winner for every fixed low-energy Yukawa vector."""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

volume = run_path(str(Path(__file__).with_name("ch39-volume-contenders.py")))


def verify(weight_modulus_error=0, gauge_prior=None):
    ctx.prec = 112
    tolerance = arb(2)**-21
    targets, _ = volume["zero_targets"](tolerance)
    x_low, x_high, logs = volume["input_bounds"]()
    top_low, top_high = logs[2].lower(), logs[2].upper()
    winner = 11737
    relative_gap = volume["q"](".0002" if weight_modulus_error or gauge_prior else ".0003")
    assert min(targets) == 10815
    if weight_modulus_error:
        assert acb.zeta_zero(min(targets)).imag > 10000
        shifted = arb(0,weight_modulus_error/10000**2)
        targets = {n:(target,weight+shifted) for n,(target,weight) in targets.items()}
        assert all(weight > 0 for _,weight in targets.values())
    nodes = 0
    for half in range(2):
        a = top_low + (top_high-top_low)*half/2
        b = top_low + (top_high-top_low)*(half+1)/2
        initial_y = logs.copy()
        initial_y[2] = a.union(b)
        data = volume["prepare"](x_low, x_high, initial_y)
        mid = initial_y[2].mid()
        stack = [(a,b,list(targets),0)]
        while stack:
            left,right,candidates,depth = stack.pop()
            nodes += 1
            top = left.union(right)-mid
            scores = {}
            for n in candidates:
                target, weight = targets[n]
                scores[n] = [(weight*lo,weight*hi) for lo,hi in
                             volume["widths"](target,data,top,tolerance,gauge_prior)]
            rivals = []
            for mode in range(2):
                lo = scores[winner][mode][0].lower()
                # A discarded rival has a certified margin on this whole slab,
                # including all of its descendants and all eight nuisance inputs.
                rivals.extend(n for n,v in scores.items()
                              if n != winner and not (1+relative_gap)*v[mode][1] < lo)
            if rivals:
                assert depth < 13, (half,left,right,rivals)
                middle = (left+right)/2
                keep = sorted(set(rivals+[winner]))
                stack.extend([(left,middle,keep,depth+1),(middle,right,keep,depth+1)])
    print("Top-coordinate boxes visited",nodes,"certified unique winner",winner)
    print("Both protocols, uniform relative advantage >", relative_gap.str(8))


if __name__ == "__main__":
    verify()

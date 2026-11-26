"""Two continuous normalizations of one identified finite height distribution."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

root=Path(__file__).resolve().parent
resolution=run_path(str(root/"ch39-phase-support-resolution.py"))
base,B=resolution["base"],resolution["B"]
parent_endpoint=run_path(str(root/"ch39-phase-support-parent.py"))["parent_endpoint"]


def main():
    weights=resolution["tables"]()[0]
    ids=sorted(weights)
    prefix=base["data"]("ch39-tick-separators.json")["prefix"]
    x={i:resolution["amb"]["upper"]["height_cell"](prefix[i-1][0]) for i in ids}
    mean=sum(x.values(),B(0))/B(len(ids))
    m=sum((weights[i]*x[i] for i in ids),B(0))
    q=sum((weights[i]*x[i]**2 for i in ids),B(0))
    # For a single distribution both scales are obtained by finite comparisons.
    assert m.hi<(q/m).lo
    for c,lo,hi in [(m,50772311235760,50772311235762),
                    (q/m,34010518992702,34010518992703)]:
        scale=c/mean
        # At each endpoint all possible positive scale factors are enclosed in
        # one pole interval before solving the secular equations.
        for mass,sign in ((lo,1),(hi,-1)):
            assert parent_endpoint(mass,ids,[weights],kinetic_scale=scale)==(sign,119)
        print("Scale",float(c.lo),"parent mass / GeV in",lo,hi,"with 119 open doublets.")
    # Mixture-affinity counterexample using two single-site records, z=1 and 3.
    assert (Q(1,2)*1**2+Q(1,2)*3**2)/(Q(1,2)*1+Q(1,2)*3)==Q(5,2)


if __name__=="__main__":
    main()

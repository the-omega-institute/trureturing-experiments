"""Minimal singlet completions, scalar degree, and determinant-phase winding."""
from itertools import combinations_with_replacement, permutations, product
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

ROOT = Path(__file__).parent


def verify():
    old = run_path(str(ROOT/"ch39-target-anomaly-completion.py"))
    indices = old["indices"]
    choices = list(combinations_with_replacement((-1,0,1),4))
    allowed = set()
    # Finite independent check of the classification throughout a charge box.
    for a,b in combinations_with_replacement(range(-5,6),2):
        charges = (a,-a,b,-b)
        for profile in choices:
            for zz in set(permutations(profile)):
                entries = tuple(sorted(zip(charges,zz)))
                if indices(entries) == (0,0,0,0,2):
                    allowed.add(entries)
    expected = {tuple(sorted(((a,1),(b,-1),(-a,0),(-b,0))))
                for a,b in product(range(-5,6),repeat=2) if a%3 == b%3 == 1}
    assert expected == allowed
    assert min(3+10*abs(a+b) for a,b in product(range(-5,6),repeat=2)
               if a%3 == b%3 == 1) == 13
    for a,b in ((1,1),(1,-2)):
        assert (-2*(a+b))%3 == 2 and (2*(a+b))%3 == 1
        print("Charged pair",a,b,"dimension",3+10*abs(a+b),
              "charged/neutral symmetric determinant windings", -2*(a+b),2*(a+b))
    ctx.prec = 150
    ext = run_path(str(ROOT/"ch39-selected-extended-flow.py"))
    v,_ = ext["initial"]()
    scale = arb("1.2209e19")
    best = v*(v/scale)**9
    old_choice = v*(v/scale)**19
    tev = v*(v/1000).root(9)
    assert arb("8.829e-25") < best < arb("8.831e-25")
    assert arb("6.385e-68") < old_choice < arb("6.388e-68")
    assert arb("1.204e16") < tev < arb("1.205e16")
    print("Unit-coefficient masses / GeV",best.str(11),old_choice.str(11),
          "dimension-13 cutoff for 1 TeV / GeV",tev.str(11))


if __name__ == "__main__":
    verify()

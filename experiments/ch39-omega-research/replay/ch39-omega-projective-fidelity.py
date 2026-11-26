"""One product-state fidelity for the shell/dyadic window, with rational Cayley coordinates."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path
import argparse

from flint import acb, arb, ctx
from sympy import Rational, expand, symbols


def circle(p, q):
    denominator = p*p + q*q
    return F(q*q-p*p, denominator), F(2*p*q, denominator)


def weight_and_fidelity(x, u):
    # x=Re(w), u=Re(D); the two radii squared are 1 and 2.
    bare = 2 * (1-x)
    squared_window = (3-2*u)*(6-4*u)*(113+72*u)/169
    weight_squared = bare**2 * squared_window
    fidelity = ((1-x)/2)**2 * ((3-2*u)/6) * ((6-4*u)/15) * ((113+72*u)/291)
    return weight_squared, fidelity


def exact():
    x, u = symbols('x u')
    score, fidelity = weight_and_fidelity(x, u)
    assert expand(419040*fidelity - 169*score) == 0
    # The highest joint Laurent coefficient certifies bidegree (2,3).
    w, z = symbols('w z', nonzero=True)
    coeff = expand(fidelity.subs({x:(w+1/w)/2, u:(z+2/z)/2})).coeff(w,2).coeff(z,3)
    assert coeff == Rational(1, 5820)
    for p, q, a, b in [(1,0,3,4),(2,17,0,1),(1,21242,4,3)]:
        x, y = circle(p,q)
        r, s = circle(a,b)
        # D=(1+i)(r+i*s) is rational complex, including at either chart's antipode.
        u, v = r-s, r+s
        assert u*u+v*v == 2 and x*x+y*y == 1
        score, kappa = weight_and_fidelity(x,u)
        assert kappa*419040 == score*169 and 0 <= kappa <= 1


def zeros(write_table):
    ctx.prec = 160
    log_two = arb(2).log()
    heights = [z.imag for z in acb.zeta_zeros(10855, 1238)]
    values = []
    grid = run_path(str(Path(__file__).with_name("ch39-omega-rational-register.py")))
    den = grid["DEN"]
    rows = ["index\tw_numerator\tdyadic_chart\tdyadic_numerator"]
    def nearest(value):
        result = (value*den+arb(1)/2).floor().unique_fmpz()
        assert result is not None
        return int(result)
    for height in heights:
        w = (acb(height, arb(1)/2) / acb(height, -arb(1)/2))
        d = arb(2).sqrt()*acb(0,height*log_two).exp()
        score, fid = weight_and_fidelity(w.real, d.real)
        direct = (abs((d-1)*(d-2)*(9+4*d)/13)/(height**2+arb(1)/4))**2
        assert direct.overlaps(score)
        assert fid.overlaps(169*direct/419040) and fid.rad() < arb('1e-39')
        k = nearest(1/(2*height))
        unit = d/acb(1,1)
        assert unit.real > 0 or unit.real < 0
        chart = int(unit.real < 0)
        slope = unit.imag/(1-unit.real) if chart else unit.imag/(1+unit.real)
        t = nearest(slope)
        rational = grid["fidelity"](k, chart, t)
        ratio = arb(rational.numerator)/(arb(rational.denominator)*fid)
        assert abs(ratio-1) < arb("2e-10")
        rows.append(f"{10855+len(values)}\t{k}\t{chart}\t{t}")
        values.append(fid)
    table = Path(__file__).with_name("ch39-omega-rational-register.tsv")
    content = "\n".join(rows)+"\n"
    if write_table:
        table.write_text(content)
    else:
        assert table.read_text() == content
    winner = 10875-10855
    assert all(values[winner] > v for j, v in enumerate(values) if j != winner)
    low, high = 10906-10855, 12054-10855
    assert all(not values[j].overlaps(values[j+1]) for j in range(low-1,high))
    assert sum(values[j-1] < values[j] > values[j+1] for j in range(low,high)) == 213
    print("Five-probe reference fidelity on 1238 arithmetic inputs; winner", values[winner].str(14))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-table", action="store_true")
    args = parser.parse_args()
    exact()
    zeros(args.write_table)

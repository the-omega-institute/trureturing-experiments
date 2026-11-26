"""Certify exact positive solutions after replacing the four phase rows by powers of the shell loop."""
from pathlib import Path
from runpy import run_path
import json
from flint import arb, arb_mat, ctx


def phase_rows():
    root = Path(__file__).resolve().parent
    phase = run_path(str(root/"ch39-shortest-triangle-clock-words.py"))
    law = run_path(str(root/"ch39-two-prime-height-law.py"))
    table = json.loads((root/"ch39-five-tick-record.json").read_text())
    den = 1 << table["bits"]
    angular_den = phase["DEN"]
    rates = law["rate_numerators"]()
    step = -43235114615559591407865398948523
    angles = []
    for _,height,_,_,_ in table["rows"]:
        n,m = (step*((v*height + (den<<23)) // (den<<24)) for v in rates)
        cells = [phase["power"](k)[0] for k in (n,m,m-n,n+m)]
        angles.append([((arb(a)/angular_den).union(arb(b)/angular_den)) for a,b in cells])
    return angles


def certify():
    ctx.prec,ctx.cap = 1664,3
    root = Path(__file__).resolve().parent
    table = json.loads((root/"ch39-five-tick-record.json").read_text())
    den = 2**table["bits"]
    angles = phase_rows()
    actual = run_path(str(root/"ch39-five-tick-record-origin.py"))["model"]
    target,g,w = actual()
    cases = [
        ("original 42 arithmetic targets",g,target),
        ("dyadic rational targets",[arb(r[1])/den for r in table["rows"]],[arb(v)/den for v in table["levels"]])]
    for label,heights,levels in cases:
        rows = [[1/(1+j*j/(h*h)) for h in heights] for j in range(42)] + list(map(list,zip(*angles)))
        mat = arb_mat(rows)
        rhs = arb_mat([[f] for f in levels] + [[0]]*4)
        center = arb_mat([[v.mid() for v in row] for row in rows])
        inverse = center.inv()
        identity = arb_mat([[int(i==j) for j in range(46)] for i in range(46)])
        error = inverse*mat-identity
        norm = max(sum(abs(error[i,j]).upper() for j in range(46)).upper() for i in range(46))
        print("residual upper",norm.str(7))
        assert norm < arb(2)**-200
        trial = (inverse*rhs).entries()
        radius = norm*max(abs(v) for v in trial)/(1-norm)
        assert min(trial) - radius > arb(".00079")
        change = max(abs(v-old) for v,old in zip(trial,w)) + radius
        print(label, "certified minimum exceeds 0.00079; change from the original weights below",change.upper().str(8))
        assert change < arb(2)**-64


if __name__ == "__main__":
    certify()

"""Generate integer matrix preconditioners from the archived finite arithmetic cells."""
import json
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx

root = Path(__file__).resolve().parent
replay = run_path(str(root/"ch39-cutoff-thresholds.py"))


def build(supports, observe=replay["observations"]):
    ctx.prec = 1000
    cases = []
    for support in supports:
        matrix, _ = observe(support)
        scale = replay["SCALE"]
        mid = arb_mat([[arb(a+b)/(2*scale) for a,b in row] for row in matrix])
        depth = 384
        entries = (mid.inv()*2**depth).entries()
        entries = [int(v.mid().floor().unique_fmpz()) for v in entries]
        n = len(support)
        cases.append({"support": support, "bits": depth,
                      "inverse": [entries[i:i+n] for i in range(0,n*n,n)]})
    return cases


if __name__ == "__main__":
    records = build(replay["supports"]())
    (root/"ch39-cutoff-threshold-inverses.json").write_text(json.dumps(records)+"\n")

"""Calibrate the input cells, or regenerate the witnesses for the integer check."""
import json
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, ctx

directory = Path(__file__).parent
output = directory/"ch39-rh-clock-integer-data.json"
counts = run_path(str(directory/"ch39-rh-clock-solution-count.py"))
fixed = counts["fixed"]


def cell(value, bits=144):
    left, right = value.lower(), value.upper()
    mantissa, exponent = right.man_exp()
    step = int(exponent)+abs(int(mantissa)).bit_length()-bits
    ulp = arb(2)**step
    return [int((left/ulp).floor().unique_fmpz())-1,
            int((right/ulp).ceil().unique_fmpz())+1, step]


def main():
    ctx.prec, ctx.cap = 1536, 3
    data = fixed["inputs"]()
    poles, rho, p, modes, k, f = data
    unit = fixed["CLOCK"]
    q02 = 1000**2*(k/22).exp()
    fields = {
        "reference": [x/unit**2 for x in poles], "rho": rho, "pivot": [p],
        "lambda": [x for x, n in modes], "q02": [q02], "target": [(15*f).exp()],
    }
    if "--write" not in sys.argv:
        stored = json.loads(output.read_text())["inputs"]
        assert set(fields) == set(stored)
        for name, values in fields.items():
            assert len(stored[name]) == len(values)
            for (lo, hi, shift), v in zip(stored[name], values):
                assert lo*arb(2)**shift < v and v < hi*arb(2)**shift
        print("All arithmetic, shell and central gauge input cells enclose their source values.")
        return
    record = {"inputs": {key: [cell(v) for v in values] for key, values in fields.items()}}
    rows = []
    for value in ("1.87758", "1.87768", "26.8343", "26.8364", "13116750", "13119360"):
        s = arb(value)
        q2 = q02/((s/8).max(arb(1))**3*(s/24).max(arb(1))).root(11)
        mf = (q2/s).sqrt()
        spectrum = fixed["shell"]["scalar_spectrum"](mf, [c/mf**2 for c in poles], rho, p, modes)
        brackets = []
        for i in range(len(modes)):
            block = spectrum[len(poles)+i*(len(poles)+1):][:len(poles)+1]
            used = 1+sum(c < q2 for c in poles)
            brackets.append([cell((x*mf**2)*(1+arb(0,"1e-17")), 96) for x, _ in block[:used]])
        rows.append({"ratio2": value, "q2": cell(q2, 96), "roots": brackets})
    for trial in (275_782_218_980, 275_782_218_981):
        sample = counts["at_unit"](data, trial)
        star = counts["at_first_pole"](sample)
        roots = counts["light_roots"](star, sample)
        rows.append({"unit": trial, "ratio2": cell((2*star).exp(), 96),
                     "roots": [[cell(x*(1+arb(0,"1e-17")), 96)] for x, n in roots]})
    record["samples"] = rows
    output.write_text(json.dumps(record, indent=1)+"\n")


if __name__ == "__main__":
    main()

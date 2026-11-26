"""Generate the integer preconditioner for the positive rational joint table."""
import json
from pathlib import Path
from runpy import run_path
from flint import arb, arb_mat, ctx


def data():
    ctx.prec=1536
    root=Path(__file__).resolve().parent
    word=run_path(str(root/"ch39-shortest-triangle-clock-words.py"))
    law=run_path(str(root/"ch39-two-prime-height-law.py"))
    record=json.loads((root/"ch39-five-tick-record.json").read_text())
    enclose=run_path(str(root/"ch39-native-positive-integers.py"))["enclose_matrix"]
    bounds=enclose(record,word,law["rate_numerators"]())
    mat=arb_mat([[arb(a+b)/(2*word["DEN"]) for a,b in row] for row in bounds])
    inv=mat.inv()
    depth=256
    nums=[[int((inv[i,j]*2**depth).floor().unique_fmpz()) for j in range(46)] for i in range(46)]
    return {"depth":depth,"numerators":nums}


if __name__ == "__main__":
    expected=json.loads(Path(__file__).with_name("ch39-native-positive-inverse.json").read_text())
    assert data()==expected


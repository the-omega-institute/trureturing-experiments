"""At most nine of the eighteen active left scalar masses can recur in the right tower."""
from fractions import Fraction
from pathlib import Path
import ast
import runpy
from types import FunctionType

from flint import arb

root = Path(__file__).resolve().parent
source = runpy.run_path(str(root/"ch39-common-source-tower.py"))
profile = source["gauge"]
q = profile["q"]
g = profile["inputs"]
model = profile["old"]
previous = profile["data"]["g"]
D,cap = model["D"],previous["cap"]
dot = previous["dot"]

def ranges(low):
    logs = [q(400).log()-x.log()/2 for x in low]
    wl = sum(logs)
    assert max(logs)<6 and wl<95
    base = [D*dot(row,g["G0"])/2 for row in previous["P"]]
    induced = previous["induced"]
    updated = dict(previous,base=base,induced=FunctionType(induced.__code__,dict(induced.__globals__,base=base)))
    beta_floor = g["beta"][0]-sum(g["D"][0])*3*cap/D
    trial = [cap,q(2),cap]
    assert g["beta"][1]+dot(g["B"][1],trial)/D < 0 and beta_floor>0
    loss = (g["d"][1]+dot(g["BD"][1],trial)/D)*wl+(8+dot(g["BF"][1],trial)/D)/3
    assert g["G0"][1]-2*loss/D>q("1/2")
    pos = previous["nonneg"]
    env = dict(g=g,m=updated,model=model,L=logs,dot=dot,cut=previous["cut"],two=q(2),four=q(4),
               Fraction=Fraction,q=q,D=D,cap=cap,
               rcR=[previous["proj"](row,model["BR"]) for row in previous["P"]],
               initial_time=pos(g["G0"][0]*D/(2*beta_floor)),
               initial_right=pos(g["G0"][0]*D/(2*q("8/15"))),
               initial_gauge=[pos(q(2)),pos(q(2)),pos(q(4))],initial_yukawa=[pos(3*cap)]*3)
    # Import just the reference inequalities, supplying the current inputs in each closure.
    path = root/"ch39-all-right-crossings.py"
    tree = ast.parse(path.read_text())
    tree.body = [n for n in tree.body if isinstance(n,ast.FunctionDef)]
    exec(compile(tree,str(path),"exec"),env)
    env["bound_color"](q("1.4"),q(4),q(2),3*cap)
    T,W = env["crossings"](logs)
    assert T>6
    smallest_ten = sum(logs[-10:])
    assert smallest_ten>W
    print("right exposure",W.lower().str(12),W.upper().str(12),"ten matching modes need",
          smallest_ten.lower().str(12))
    return T,W

if __name__=="__main__":
    ranges(profile["low"])
    ranges(source["mass_bounds"]())

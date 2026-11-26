"""Lower genus of the physical compactified mass curve from disjoint ramification clusters."""
from pathlib import Path
from runpy import run_path
from flint import acb, arb

branch = run_path(str(Path(__file__).with_name("ch39-shell-gw-mass-transposition.py")))
q = arb
centers = {
    43: [("-4.2012188","-1.1576734"),("-3.5282969","2.6033374"),
         ("-.1262802","-4.7384616"),("-3.4176287","-3.5569844"),
         ("-3.7080711","3.5245817"),("-.0315593","5.2806656"),
         ("-3.2948645","-4.6988058")],
    195: [("-4.2012288","-1.1576201"),("-3.5283475","2.6033448"),
          ("-.1262131","-4.7385144"),("-3.4176839","-3.5570217"),
          ("-3.7081359","3.5246154"),("-.0313690","5.2807122"),
          ("-3.2949131","-4.6988816")]
}


def prove(index):
    data = branch["symbol"](index)
    trial = [acb(x,y) for x,y in centers[index]]
    for i,c in enumerate(trial):
        assert abs(c)>4 and all(abs(c-z)>q(".1") for z in trial[:i])
        branch["certify"](index,c,q(".00004"),data)
    # Each non-unit-circle cluster is mirrored by the physical real structure.
    ramification = 2*len(trial)*52
    lower = 1+(ramification-2*53)//2
    print("certified disjoint simple branch values",ramification,"genus at least",lower)


if __name__=="__main__":
    for n in (43,195):
        prove(n)

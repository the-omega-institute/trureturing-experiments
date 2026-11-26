"""Determine the minimum common positive integer tick for the first 54 normalized RH readings."""
from pathlib import Path
from runpy import run_path
from flint import arb
import json


def main():
    root=Path(__file__).resolve().parent
    source=run_path(str(root/"ch39-native-five-tick-boundary.py"))
    certify=run_path(str(root/"ch39-prime-clock-height-horizon.py"))["main"]
    data=json.loads((root/"ch39-native-eight-tick.json").read_text())
    for tick in (1,2,3,4,6,7):
        certify(tick,15)
        c=data[str(tick)]
        # For a box bound, arbitrary independent cosines are allowed in the high tail.
        # For ticks one and seven, use the exact two-angle dependence there.
        peak={1:".369381",7:".925928"}.get(tick)
        source["separate"](c,tick,44,".003",".002",peak,96 if tick==1 else 272)
    source["main"]()
    certify(8,15)
    support=data["8"]["support"]
    assert len(support)==58
    run_path(str(root/"ch39-native-clock-held-out.py"))["main"]([(support,54)],arb(".00037"),8)
    print("The least common positive integer multiplier at 54 consecutive readings is eight.")


if __name__=="__main__":
    main()

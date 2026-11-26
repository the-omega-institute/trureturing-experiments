"""Exact finite separation and positive integer votes for 45 time offsets and four phases."""
import json
from pathlib import Path
from runpy import run_path

root = Path(__file__).resolve().parent


def main():
    load = lambda name: json.loads((root/name).read_text())
    record = load("ch39-adjacent-offset-certificate.json")
    case = record["case"]
    assert case["offsets"] == [*range(43),1682,1683] and case["cut"] == 3492
    assert len(set(case["support"])) == len(case["support"]) == 49
    prefix = load("ch39-tick-separators.json")["prefix"]+record["tail"]
    lows = load("ch39-minimum-tick-finite-inputs.json")["levels"]
    values = lows[:43]+[load("ch39-fixed-register-extra-target.json")["floor"],record["target1683"]]
    assert len(prefix)==3492 and prefix[-1][0]+1<2**415
    assert all(x[0]+1<y[0] for x,y in zip(prefix,prefix[1:]))
    run_path(str(root/"ch39-prime-clock-height-horizon.py"))["main"](5,15)
    verify = run_path(str(root/"ch39-extra-offset-cutoffs.py"))
    verify["replay_lower"](case,values,prefix)
    verify["replay_upper"](case,values,prefix)
    print("Five ticks, 45 time readings and four phases: least largest zero index 3492, with 49 positive heights.")


if __name__ == "__main__":
    main()

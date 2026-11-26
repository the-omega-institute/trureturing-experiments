"""Least highest zero after the next ordinary held-out RH time is enforced."""
import json
from pathlib import Path
from runpy import run_path
import sys

root=Path(__file__).resolve().parent
path=root/"ch39-fifty-site-repair.json"


def main():
    if len(sys.argv)==3 and sys.argv[1]=="--write":
        case=json.loads(Path(sys.argv[2]).read_text())
        a,case["weights"]=run_path(str(root/"ch39-extra-offset-cutoffs-origin.py"))["compute"](case)
        assert a==case["a"]
        path.write_text(json.dumps(case)+"\n")
    case=json.loads(path.read_text())
    if "--origin" in sys.argv:
        a,w=run_path(str(root/"ch39-extra-offset-cutoffs-origin.py"))["compute"](case)
        assert a==case["a"] and w==case["weights"]
    assert case["cut"]==2419 and case["offsets"]==[*range(44),1682]
    assert len(set(case["support"]))==len(case["support"])==49
    def data(name):
        return json.loads((root/name).read_text())
    prefix=data("ch39-tick-separators.json")["prefix"]
    y=data("ch39-minimum-tick-finite-inputs.json")["levels"][:44]+[data("ch39-fixed-register-extra-target.json")["floor"]]
    assert prefix[2418][0]+1<2**415
    run_path(str(root/"ch39-prime-clock-height-horizon.py"))["main"](5,15)
    checks=run_path(str(root/"ch39-extra-offset-cutoffs.py"))
    checks["replay_lower"](case,y,prefix)
    checks["replay_upper"](case,y,prefix)
    print("Including j=43: least highest index 2419, attainable at 49 positive heights.")


if __name__=="__main__":
    main()

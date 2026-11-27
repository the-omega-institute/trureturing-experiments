"""第八壳色阈值二的可行参数下界。"""
from pathlib import Path
from runpy import run_path

from flint import arb


CENTERS = (
    "75.16623387842753",
    "7.472616033279199",
    "77.0924672089881",
    "64.41308758794867",
    "7.19432327017437",
    "7.471975487008222",
    "0.24443833381882607",
    "40.368451868210826",
    "7.480226948737031",
    "0.24443833381882607",
    "77.90477970362736",
    "36.010486570983815",
    "0.24443833381882607",
    "41.92910075130166",
    "7.47229016202135",
    "0.24443833381882607",
    "55.94487695916509",
    "34.89714491069651",
    "7.47363894762853",
)
ACTIVE = (15, 10, 16, 8, 7)


def verify(root):
    case = run_path(str(root / "ch39-eighth-arc-positive-probe.py"))
    t = arb(".736102") + arb(0, arb("1e-12").upper())
    case["verify"](root, parameter_case=(t, CENTERS, ACTIVE))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)

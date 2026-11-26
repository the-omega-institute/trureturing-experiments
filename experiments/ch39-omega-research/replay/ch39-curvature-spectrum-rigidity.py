"""Constant scalar curvature couplings with the original charged and neutral spectra."""
from pathlib import Path
import runpy
from flint import arb

shell = runpy.run_path(str(Path(__file__).with_name("ch39-short-arc-vacuum.py")))
q = shell["q"]
r = shell["r"]
ranks = shell["count"]
# Higgs shifts are given once per vertex orbit in M^2 units, and the
# right-hand internal connection has rank two. The scalar graph has rank one.
higgs_trace = 2 * sum(n * (v + p) for n, v, p in
                      zip(ranks, shell["potential"], shell["shiftbox"])) / r**2
assert len(ranks) == len(shell["potential"]) == len(shell["shiftbox"]) == 5
assert sum(ranks) == 19
assert q("120.05") < higgs_trace < q("127.95")
weyl = r * sum(shell["total"]) / (2 * arb.pi().sqrt())
assert q("78.44811") < weyl < q("78.44813")

if __name__ == "__main__":
    print("neutral rank-one Weyl coefficient", weyl.str(14))
    print("trace of the charged Higgs potential in Q units",
          higgs_trace.lower().str(14), higgs_trace.upper().str(14))

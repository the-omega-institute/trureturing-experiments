"""Reconstruct the exact certificate and retain its mathematical arrays, not logs."""
from hashlib import sha256
from pathlib import Path
import json
import runpy

import numpy as np

HERE = Path(__file__).resolve().parent
EXPERIMENTS = HERE.parent
state = runpy.run_path(str(HERE / "independent_precheck.py"))
arrays = {
    "matrix.u16le": state["m"],
    "pivot-rows.u16le": state["pivots"],
    "denominators.u16le": state["denominators"],
    "input-wedges.u16le": state["inputs"],
}
data = HERE / "data"
data.mkdir(exist_ok=True)
records = {}
for name, values in arrays.items():
    array = np.asarray(values, dtype="<u2")
    payload = array.tobytes(order="C")
    (data / name).write_bytes(payload)
    records[name] = {"shape": list(array.shape), "sha256": sha256(payload).hexdigest(),
                     "bytes": len(payload)}
sources = [
    "auric-shared-d-frame-central-pb/frame_kernel.py",
    "auric-shared-d-frame-central-pb/central_pb.py",
    "auric-shared-d-frame-central-pb-precheck/independent_precheck.py",
    "auric-shared-d-frame-central-pb-precheck/retain_certificate.py",
]
result = {
    "scalar_field": "Q (real and imaginary parts of Q(i) before reduction)",
    "prime": state["P"], "i_image": 469,
    "array_encoding": "unsigned 16-bit little-endian, C row-major, zero-based indices",
    "arrays": records,
    "source_sha256": {name: sha256((EXPERIMENTS / name).read_bytes()).hexdigest()
                      for name in sources},
    "rank_mod_1009": state["rank"],
    "maximal_minor_mod_1009": state["determinant"],
    "minor_rows": "leftmost pivot columns of RREF(matrix transpose), increasing order",
    "minor_columns": "all columns in original order",
    "input_wedge_rank_mod_1009": 105,
    "full_map_dimension": 1680, "central_kernel_dimension": 105,
    "householder_denominator_count": len(state["denominators"]),
    "householder_denominator_minmax": [min(state["denominators"]), max(state["denominators"])],
    "central_witnesses": {
        "a=1": {"B_plus=B_minus": "(80*s^6+12*s^4+6*s^2+1)/576",
                "L": "-(44*s^2+1)/12",
                "L^2-4*B": "s^2*(41+962*s^2-40*s^4)/72",
                "small_range": "0<abs(s)<=1/4"},
        "a=2": {"B_plus=B_minus": "(2*s^6+6*s^4+9*s^2+1)/9",
                "L": "-2*(5*s^2+1)/3",
                "L^2-4*B": "4*s^2*(1+19*s^2-2*s^4)/9",
                "small_range": "0<abs(s)<=1"},
        "large_range": "abs(s)>a^2/4: x=0, L=-2*s, B_plus=B_minus=a^4/16",
        "cases": "a=1,2; sigma=-1,+1; phase=1,i; both real-Gram budgets",
    },
}
(data / "certificate.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
print(json.dumps({"retained_arrays": records, "maximal_minor_mod_1009": state["determinant"]}, sort_keys=True))

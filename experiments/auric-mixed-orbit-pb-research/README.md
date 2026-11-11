# Mixed-orbit calibrated raw PB checks

With Python 3.12 or newer, install `requirements.txt` and run:

```sh
python3 experiments/auric-mixed-orbit-pb-research/certificate.py
```

`certificate.py` regenerates `evidence.json`: exact symbolic exterior/Gram and
weighted Cauchy–Schwarz identities, the orbit-spanning determinant, the alpha=2
pure failure, and deterministic finite numerical diagnostics. Its arbitrary-real
lambda identities support the analytic scalar-family proof in the main theory
volume `AURIC_FULL_PB_RAW_CALIBRATION_AND_COMMON_SOURCE_TRANSPORT_CONTINUATION.md`
§52. Finite diagnostics do not prove the universal alpha/U/V/q/x claim.

The unchanged calibrated definitions are the three-actor kernel-lift volume
§41.1, source pin `17245b37e7cf1381419fb570b94b23b421efd1fd`. The actor label is
h=K_sigma(c,x); only g=D_alpha(c,t) is varied. This is a matrix-tests-only
obstruction to replacing positive-cone tests by signed tomography, with no
native Omega countermodel or complete quantum-origin claim.

Standard informational completeness of a nontrivial unitary orbit is reused
from D'Ariano–Perinotti–Sacchi, *Informationally complete measurements and groups
representation*, arXiv:quant-ph/0310013, §4.2, equations (29)–(31).
Only program source and result data are stored here; the full PB argument and
its source boundaries are in the owning main volume.

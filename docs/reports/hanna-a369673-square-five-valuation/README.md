# Hanna A369673: square classification and five-adic orders

`check.py` uses Python 3 integer arithmetic to check the full four-clause
conjecture package in [Paul D. Hanna's OEIS A369673](https://oeis.org/A369673),
revision #21, together with the auxiliary factorial-valuation definition in
[OEIS A127428](https://oeis.org/A127428), revision #21. The source fidelity
checks reproduce all 15 displayed A369673 terms and the first 31 displayed
A127428 terms directly from their defining expressions.

Run `python3 check.py`; the deterministic output is in `output.txt`.
The conjectures are checked for every product index through 400: C.1 for
N = 0..400, C.2 for n = 0..199, C.3 for n = 0..100, and C.4 for n = 0..99.
The program also checks the two product pairings, the exact lifting identity
for j = 1..5000, and the odd-factor factorial bridge for r = 0..201.

All arithmetic, square roots, factorials and valuations are exact; no
floating point or external package is used. A nonpositive valuation argument
is rejected. Zero-index and empty-product conventions are included. These
are finite checks supporting a preregistration, and do not prove a universal
statement or supply a Lean proof.

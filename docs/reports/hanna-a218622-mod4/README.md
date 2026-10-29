# Hanna A183161/A218622: residues modulo four

Run `python3 check.py` with Python 3 and SymPy installed. The check uses exact
integers, rational fractions, and polynomial arithmetic over finite rings.
The stored output is `output.txt`; success ends with `ALL_CHECKS_OK`.

The literal A183161 finite sum, its normalized self-convolution recurrence,
and all 301 published integer coefficients agree at indices 0 through 300.
The proposed 15-state automaton agrees with all 2049 published A218622 residues.
All 30 proposed numerator-pair transitions are verified by both sparse modular
multiplication and independent dense integer convolution. The resultant and
rational-diagonal numerator identities are checked exactly. Further finite
checks compare the generating-function bridge through degree 150 and directly
divide the two bivariate rational series through diagonal index 64.

These are finite source-fidelity and certificate checks. The formal power-series
bridges, section theorem and binary induction have no Lean proof here. A compiled
statement definition is not a proof of the conjecture.

Sources: Paul D. Hanna, OEIS [A183161](https://oeis.org/A183161), revision 31,
and [A218622](https://oeis.org/A218622), revision 9. Their conjecture is dated
November 3, 2012. The unmodified source tables are
[b183161.txt](https://oeis.org/A183161/b183161.txt) and
[b218622.txt](https://oeis.org/A218622/b218622.txt), attributed in OEIS to
Paul D. Hanna. These OEIS data are licensed under
[CC BY-SA 4.0](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement).
General methods: Rowland–Yassawi, *Automatic congruences for diagonals of rational
functions*, arXiv:1310.8635v2, DOI 10.5802/jtnb.901; the specific polynomial
certificate is the candidate supplied for independent checking.

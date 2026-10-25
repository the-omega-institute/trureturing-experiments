# Hanna A389472, A389474 and A389475 congruences

Run `python3 check.py`. Python 3's arbitrary-precision integers suffice; there are no dependencies or network requests. `output.txt` is the complete standard output; a successful run exits zero and ends in `ALL_CHECKS_OK`.

The source is Paul D. Hanna's OEIS entries [A389474](https://oeis.org/A389474) (revision 11), [A389475](https://oeis.org/A389475) (revision 8), and [A389472](https://oeis.org/A389472) (revision 7). Their displayed DATA terms are embedded as fidelity fixtures. The script independently constructs the integer coefficients from the defining equations, checks each division by two is exact, and reproduces all 24, 23 and 41 displayed terms, respectively.

A389474 is constructed by the triangular coefficient recurrence for `A(x²+2x³+2A²)=3A²`; its composition equation is independently recomputed by truncated Horner evaluation. A389475 is constructed by the triangular recurrence for `A(3x²)=2x²+A²−2A³`, then checked directly. Both series are generated through coefficient 400. The script tests A389474 (3), A389474 (4), A389475 (4.a), and A389475 (4.b) at every coefficient through degree 400. For A389474 (4), all relevant powers `k=1..134` are checked, and larger powers vanish in this truncation because each base has order three.

A389472 uses the binomial recurrence with seeds `a(0)=0`, `a(1)=a(2)=1`, matching the source program and the existing repository owner. It is constructed through coefficient 1202; the defining equation is checked independently by binomial composition through degree 1202. Its mod-3 conjecture is checked for the 400 target indices `n=2..401`. The parity conjecture, already settled in the main repository, is also a regression check and is excluded from new preregistration.

These are exact finite checks, not proofs of the universally quantified conjectures. The output records SHA-256 fingerprints of the generated coefficient arrays using compact JSON, including the constant coefficient. Source numeric evaluations at nonzero real arguments are outside this check.

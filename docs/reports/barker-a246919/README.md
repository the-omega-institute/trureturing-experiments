# Barker's maximal integral cevian generating function

The exact checks support the generating-function conjecture in [OEIS A246919](https://oeis.org/A246919), FORMULA, attributed to Colin Barker (6 June 2016). The source snapshot is revision 27 (28 June 2026). A cevian may meet the extended opposite side: its coordinate x is any integer except 0 and n. Its positive integral length y satisfies x² − y² − nx + n² = 0.

Run with Python 3 and SymPy:

```sh
python3 docs/reports/barker-a246919/check.py
curl -L https://oeis.org/A246919/b246919.txt -o /tmp/b246919.txt
python3 docs/reports/barker-a246919/check.py --bfile /tmp/b246919.txt
```

The stored output is from the second invocation, with Python 3.14.8 and SymPy 1.14.0. The b-file snapshot's SHA-256 is 7d7f15f0a8a94f25ecd29188cc8e0860793cf5cc9977049167e9b592440ab6e6.

The check reproduces all 53 printed terms by an integer translation of the source's PARI program. It independently enumerates every admissible factor pair for n = 1,…,300, cross-checks complete solution sets by scanning the original quadratic for n = 1,…,80, and compares exact rational-function coefficients with the proposed maximum through n = 10000. Every one of the 10000 b-file values agrees. A separate exact symbolic calculation reduces the difference of the residue-class generating functions and the printed rational function to zero.

The proposed maximum is zero at n = 1, 2, 4; otherwise it is (3n²+1)/4 for odd n, (3(n/2)²+1)/2 for n ≡ 2 modulo 4, and 3(n/4)²+1 for n divisible by 4. The unbounded argument still requires proving these maxima from the Diophantine definition. These finite checks are corroboration and do not constitute a Lean proof.

A separate source error is recorded: the printed recurrence claims n > 14, but its right side at n = 16 is 45, while the exact maximum is 49 with witness x = 55. This does not refute the generating-function conjecture. The proposed correction is n ≥ 17.

Initial sequence values and the PARI algorithm are attributed to Colin Barker and the OEIS contributors; see the [OEIS End-User License Agreement](https://oeis.org/LICENSE). The check implementation is supplied under this repository's license.

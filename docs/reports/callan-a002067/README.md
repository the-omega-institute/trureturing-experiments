# Exact checks for Callan's A002067 conjectures

Source: [OEIS A002067](https://oeis.org/A002067), revision 116, two COMMENTS conjectures of David Callan (July 21, 2017). The recurrence has initial value a(0) = 1. Positions and vertex labels start at one, and an upstep uses the height of its upper endpoint. The empty path and empty matching each have weight one.

Run with Python 3 (standard library only):

```sh
python3 docs/reports/callan-a002067/check.py
```

`output.txt` contains the output, ending in `ALL_CHECKS_OK`. The checks use integers and rational numbers throughout.

- Reproduce the source example UUDUUDDD, with positions 1, 2, 4, 5 and upper heights 1, 2, 2, 3, giving weight 480. Reproduce the first eleven terms and Dominici's first three derivative polynomials.
- Independently calculate the NAME recurrence, the perfect-matching weighted sum, every Dyck path with the source weight, and P(2n, 0), for 0 <= n <= 10.
- Explicitly enumerate every matching for 0 <= n <= 6; verify every matching-to-path fibre has the product of upstep upper heights as its cardinality and the stated total weight.
- Enumerate every nonnegative walk of each length through 12 and compare every coefficient with the polynomial recurrence.
- Use rational formal-series arithmetic through degree 20 to check y' = exp(y^2/2), A exp(-y^2/2) = 1, and C'' = A through degree 18, where C = 1 - exp(-y^2/2).

The eleven common values are 1, 1, 7, 127, 4369, 243649, 20036983, 2280356863, 343141433761, 65967241200001, 15773461423793767. These are finite checks supporting a proposed general argument; they do not prove either universal conjecture. No Lean proof is included.

Relevant published ingredients are Callan, *A combinatorial survey of identities for the double factorial*, arXiv:0906.1317, Section 3, Proposition 1, and Dominici, *Asymptotic analysis of the derivatives of the inverse error function*, arXiv:math/0607230v2, Section 3. The check program implements the stated recurrences and elementary enumeration independently; it contains no copied third-party program.

# Ceiling Pisot L(6,7): exact finite checks

`python3 check.py` requires Python 3 and its standard library. It prints
`ALL_CHECKS_OK` and exits 0 when every check passes. `output.txt` records the
result for 10,002 terms. `--terms N` changes the finite range.

The program independently generates OEIS A048585 by rational ceiling and by
`v(0)=3, v(1)=4, v(2)=6, v(n+3)=v(n+2)+v(n)`, followed by `U=v+3`.
It reproduces every displayed source term and checks the proposed order-four
identity, the second-difference and determinant recursions, their initial
values, and the ceiling error window. All evidence uses integers and
`fractions.Fraction`; no floating point is used.

Source: David W. Wilson, [OEIS A048585](https://oeis.org/A048585/internal),
revision 18. The displayed sequence data are attributed to Wilson and the
OEIS Foundation under [CC BY-SA 4.0](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement).
The independent check program follows this repository's license.

These computations certify only the stated finite ranges. They do not prove
the all-index assertion, certify literature novelty, or constitute a Lean proof.

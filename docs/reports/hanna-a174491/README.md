# Hanna A174491 denominator checks

For the formal power series

\[
F(x)=\exp\left(\sum_{m\ge1}\frac{x^m}{m2^{m^2}}\right)
     =\sum_{n\ge0}c_nx^n,
\]

OEIS A174491 conjectures that the reduced denominator of \(c_n\), for every
positive integer \(n\), is \(2^{n^2+v_2(n)}\). The constant coefficient is 1.

Sources: [A174491](https://oeis.org/A174491), revision 7, and
[A174490](https://oeis.org/A174490), revision 7. The entries are authored by
Paul D. Hanna and explicitly mark the denominator formula as a conjecture.

Run `python3 check.py` with Python 3.10 or newer. The program uses only the
standard library, integers, and exact rational arithmetic. Its output is in
`output.txt`.

The checks reproduce all published denominators through index 11 and all
published numerators through index 12. A coefficient recurrence and a product
of exponential factors agree through index 24. The proposed denominator
identity holds through index 80. Partition terms through index 24 have their
unique smallest 2-adic valuation at the one-part partition. Odd-prime
logarithmic differences and the required lifting congruences are checked in
the finite ranges printed by the program.

These computations provide source fidelity and finite checks of a proposed
general argument. They do not prove the conjecture for all indices.

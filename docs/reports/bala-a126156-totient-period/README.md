# A126156: mod-4 obstruction to Bala's eventual totient-period conjecture

`python3 check.py` uses only Python 3 standard-library integers and `fractions.Fraction`. Its expected output is `output.txt`, ending in `ALL_CHECKS_OK`.

The program constructs the sequence independently from the A126155 integer triangle, the tangent differential recurrence, and the rational formal-series identity `F(x)^2*cos(sqrt(2)*x)=1`, with `F(0)=1`. All three agree for indices 0 through 200. Fidelity checks reproduce all 17 displayed A126156 terms, the A126155 row at index 4, and Bala's 28 displayed residues modulo 55. Shift 20 modulo 55 is checked at indices 1 through 180.

The mod-4 checks verify the row invariant, oddness, the recurrence `a(n+1) = (2n+1)*a(n) mod 4`, and the persistent two-step difference `a(n+2)-a(n) = 2 mod 4` on the stated finite ranges. This gives the cycle `1,1,3,3` in every tested block, while `phi(4)=2`. Finite computation supports the separately preregistered induction argument; it does not exclude every onset by itself and is not a Lean proof.

Sources: Peter Bala's COMMENTS conjecture in [OEIS A126156](https://oeis.org/A126156), revision 110, and its defining e.g.f. attributed to David Callan; the integer-triangle algorithm in [OEIS A126155](https://oeis.org/A126155), revision 26. Sequence terms and the mod-55 example are attributed to those entries, whose content is available under [CC BY-SA 4.0](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement). The mathematical statement has no source arXiv identifier.

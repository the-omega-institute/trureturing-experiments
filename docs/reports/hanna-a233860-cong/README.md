# Hanna A233860: exact congruence checks

Run `python3 check.py`; the deterministic output is in `output.txt`.
Python 3 standard-library integers and `fractions.Fraction` supply all arithmetic.
The program reproduces all 20 DATA terms of [OEIS A233860](https://oeis.org/A233860)
(revision 10), computes 160 exponential coefficients with the integer derivation
`D = P(t) d/dt`, and checks every applicable parity, 5-power and 7-power
congruence in that range. An independent rational reversion through degree 32
checks the defining ODE and both inverse compositions. The degree-24 inverse
certificate over F_5 and all ten exceptional low-degree coefficients are checked
exactly. The output is finite evidence; it does not certify an unbounded theorem.

Sources: Paul D. Hanna's [internal entry](https://oeis.org/A233860/internal),
its [JSON](https://oeis.org/search?q=id:A233860&fmt=json), and Vaclav Kotesovec's
[table](https://oeis.org/A233860/b233860.txt). The program uses the entry's DATA
and the ten low-degree table values as fidelity references. The check program
is independently written and uses no third-party Python code.

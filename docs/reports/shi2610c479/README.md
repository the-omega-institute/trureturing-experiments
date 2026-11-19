# Shi–Yau–Zuo Conjecture 4.79

Run `python3 check.py` with Python 3.12 or later. The standard library suffices. The output uses exact integer arithmetic and modular polynomial coefficients.

The first check constructs the 16-element square-zero algebra over F₂. It independently generates ideals by closing their generators under addition and multiplication by ring elements. The principal-ideal map and the proposed three-component map reproduce Example 4.5(1) and Proposition 4.67 of the source on this algebra. The proposed map gives a T-full maximal ideal absent from its list of value ideals.

The second check computes the dependence witness and its factorization over the prime fields of orders 2, 3, 5 and 7. The nonzero polynomial is divisible by every coordinate factor, and each rational pair kills at least one factor. These are finite tests of the family construction.

[Shi, Yau and Zuo, Izvestiya Mathematics 88 (2024), 1185–1220](https://doi.org/10.4213/im9488e).

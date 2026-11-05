# Even normalized Eulerian maxima log-convexity

`check.py` uses only integer arithmetic and exact `fractions.Fraction` values. It independently enumerates permutations for rows 1 through 8, cross-checks the Eulerian recurrence, verifies the normalized maxima `f(4)=11/24`, `f(6)=151/360`, and `f(8)=15619/40320`, and checks all centered inequalities through `p=100`. It also checks the ratio bounds and polynomial separation used in the preregistration outline.

The v3 TeX source has a table of initial `f(m)` values that is inconsistent with its own permutation definition (for example it prints `f(4)=23/48`, whereas the defined maximum is `A(4,3)=11`, hence `f(4)=11/24`). The table is not used as evidence; the definition, exact enumeration, and recurrence are used instead. The journal version was not accessible from the publisher during this screen.

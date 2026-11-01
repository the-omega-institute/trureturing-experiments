# OEIS A236401: Kurepa determinant asymptotic

`check.py` uses SymPy integer determinants and Python exact integers and fractions. It reproduces every published OEIS value for n from 7 through 25, verifies the derangement recurrence and its finite alternating-sum expression, and checks the proposed all-index identity for the determinant through n equal to 40:

K_n = (-1)^(n+1) (6 (n-5)! + d_(n-4) + d_(n-5)).

The output is stored in `output.txt`. These are finite exact checks; they do not replace the all-n proof or its limiting argument.

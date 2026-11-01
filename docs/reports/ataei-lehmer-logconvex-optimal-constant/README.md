# Ataei log-convex optimal constant

`check.py` uses only Python `fractions.Fraction` arithmetic. It implements the source's piecewise-affine primitive
\(F(s)=n s-n(n+1)/2\) for \(n=\lfloor s\rfloor\), then checks the proposed extremal pair
\(\varphi_2=\eta F(s)\) and \(\varphi_1=\eta(F(s+1/2)-s/2)\).

For several exact positive values of \(\eta\), it verifies the common recurrence increment, the exact growth ratio, the periodic difference
\(p(s)=-(\eta/2)\min\{\{s\},1-\{s\}\}\), and the lower bound \(\eta/8\) after any constant normalization. The final line is `ALL_CHECKS_OK`.

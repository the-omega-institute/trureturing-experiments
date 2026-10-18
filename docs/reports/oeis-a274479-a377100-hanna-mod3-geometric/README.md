# OEIS A274479 and A377100: a(n) = 1 (mod 3)

`check.py` solves both functional equations over Q to order N (default 200; the reading below uses 160) and checks the settlement steps of the open-problem preregistrations in `the-omega-institute/trureturing`:

- **A274479** (A(x)^2 = A(x^2/(1-2x-4x^2))). The reversion R satisfies R(y^2)(1 - 2R - 4R^2) = R^2 and its recursion divides only by 2, exactly. Modulo 3, x/(1-x) solves the equation.
- **A377100** (A(x) = A(x^3)/A(x^2) + A(x)^2). The recursion is division-free. Modulo 3, x/(1-x) solves the equation.

Run: `python3 check.py 160` (exit 0, `RESULT PASS`; output in `out.txt`).

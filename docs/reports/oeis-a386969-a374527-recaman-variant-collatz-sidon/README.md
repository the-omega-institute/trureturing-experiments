# OEIS A386969 and A374527: a Recamán variant and Collatz trajectories that are Sidon sets

**A386969.** This is the Recamán variant a(n) = a(n−1) − n if nonnegative and new, else a(n−1) + n + 3. Two equivalent conjectures describe it:
- Z. He's generating function x(4+2x+11x³−2x⁴−9x⁵)/(1−x²)²;
- B. McEachen's closed forms a(n) = n + 5 (odd n > 2) and 2n + 7 (even n > 3).

They follow by induction. For n ≥ 4 the earlier values are {0, 2, 4} ∪ {m+5 : m odd} ∪ {2m+7 : m even}:
- at even n the backward step lands on 4, which is already used;
- at odd n it lands on the new even value n + 5.

**A374527.** This lists the odd numbers whose Collatz trajectory is a Sidon set. R. Stephan conjectures that it contains every (4^k − 1)/3 except 5. The trajectory of x = (4^k − 1)/3 is x followed by 4^k, 2^(2k−1), …, 1. For k ≥ 3 the sums x + 2^a and 2x have binary weight ≥ 3, while sums of two powers of two have weight ≤ 2. For k = 2, 5 + 5 = 8 + 2.

`check.py` (Python 3) checks:
- the published data;
- the closed forms, the g.f. identity (1−x²)²A(x) = numerator, and the set invariant for n ≤ N;
- the Sidon property and the trajectory shape for k ≤ K.

Run `python3 check.py 3000 60`. Exit 0 and `RESULT PASS` mean every check holds.

# OEIS A355257, A355407, A355414: Tebni's congruences for log((1−x)/(1−2x))/(1−x)^n

A355257 is the array A(n, k) = k!·[x^k] log((1−x)/(1−2x))/(1−x)^n. A355407 and A355414 are its rows n = 4 and n = 5. The conjectures checked are:
- **C1.** A(n, p) ≡ −1 (mod p) for p prime. This is stated in A355257, A355407 and A355414.
- **C2.** For k ≥ 1, k ≠ 4: k | A(n, k) iff k is not prime.
- **C4.** For p prime and n even, A(n, p)/(p − 1) is an integer ≡ 1 (mod p).

All three follow from A(n, k) = Σ_{j=1..k} (k!/j)(2^j − 1) C(k−j+n−1, n−1):
- the j < k terms are divisible by k;
- the j = k term is (k−1)!(2^k − 1), which handles C1 by Wilson and Fermat, and C2 because composite k ≠ 4 divides (k−1)!;
- modulo p − 1 only the j = p − 1 term survives, which handles C4.

`check.py` (Python 3 with sympy) does three things:
- it compares the formula with a series expansion;
- it checks the published data of all three entries;
- it checks C1, C2 and C4, and the reduction modulo p − 1.

Run `python3 check.py 30 200`. Exit 0 and `RESULT PASS` mean every check holds.

# OEIS A384531 and A367865: two multiplicative gcd-sum identities

Two conjectured identities between multiplicative functions:
- **A384531** (W. Schulte). a(n) is multiplicative with a(p^e) = ((2e+1)p − 2e)·p^(e−1), and the conjecture is a(n) = Σ_{i=1..n} gcd(i, n)·gcd(i+1, n).
- **A367865** (W. I. Hurt; conjecture by V. Yanev). a(n) = Σ_{d|n} d·φ(d)·μ(d)², and the conjecture is a(n) = (rad(n)/n)·Σ_{k=1..n} rad(nk)/rad(k).

`check.py` (Python 3 with sympy) checks:
- both conjectures exactly, the second in rational arithmetic;
- the published data;
- the two local facts the proofs use: at most one of gcd(i, p^e) and gcd(i+1, p^e) exceeds 1, and rad(nk)/rad(k) = rad(n)/rad(gcd(n, k)).

Run it with `python3 check.py 1500`. Exit 0 and `RESULT PASS` mean every check holds. `run.log` is the output.

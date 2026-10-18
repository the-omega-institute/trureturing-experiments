# OEIS A377106 (Hanna): a(n) is odd iff n is a power of 2

Entry: "G.f. A(x) satisfies A(x)^3 = A( x^3 + 9*x*A(x)^3 ), with A(0)=0, A'(0)=1." Conjecture: "a(n) is odd iff n = 2^k for some k >= 0."

Route: the series reversion R = A⁻¹ satisfies R(y³) = R(y)³ + 9y³R(y).
- **Integrality.** The coefficient recursion is 3r_n = [y^{n+2}](R(y³) − R_{<n}(y)³ − 9y³R_{<n}(y)), and the bracket is divisible by 3 because R_{<n}(y)³ ≡ R_{<n}(y³) (mod 3). So R, and hence A, has integer coefficients.
- **Parity.** R ≡ y + y² (mod 2), and the compositional inverse of y + y² over GF(2) is Σ_k x^{2^k}.

## Program

`check.py` (standard library only). Run `python3 check.py`. It works to order 150 and prints `RESULT PASS` when every check holds:
- every division by 3 in the recursion for R is exact;
- R ≡ y + y² (mod 2);
- the inverse A agrees with the first 12 OEIS terms;
- A satisfies A³ = A(x³ + 9xA³);
- a(n) is odd exactly for n ∈ {1, 2, 4, …, 128}.

Reading: `run.log` (exit 0, `RESULT PASS`, 1.7 s).

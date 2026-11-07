# Yin finite-polygamma majorization checks

Source: Li Yin, “Monotonic properties for ratio of the generalized (p,k)-polygamma functions”, Journal of Mathematical Inequalities 16(3) (2022), 915–921, DOI 10.7153/jmi-2022-16-62, Conjecture 4.1. The publisher PDF is the source; there is no arXiv identifier.

Run `python3 docs/reports/yin-polygamma-majorization/check.py` with Python 3 and SymPy 1.14.0. The program uses symbolic algebra and `fractions.Fraction`, with no floating-point evidence. The saved standard output is `check-output.txt` and ends in `ALL_CHECKS_OK`.

Fidelity checks differentiate the source’s finite gamma expression (1.4) to recover the finite digamma expression (1.5). Further checks differentiate (1.5), explicitly exposing the incompatible finite-sum sign printed in (1.6). The derivative definition and (1.5) govern the ratio.

For p = k = 1, alpha = (2,2,2), beta = (1,2,3), the exact quotient increases from 54/85 at x = 1 to 2450/3783 at x = 2, a difference of 3968/321555. The vectors satisfy strict majorization and every derivative order is positive.

The family argument uses S_r(x) = sum from j = 0 to p of (x+jk)^(-r). For p >= 1, k > 0, m >= 2, the quotient is R = [m/(m+1)] S_(m+1)^2/(S_m S_(m+2)). Its strict deficit follows from

    S_m S_(m+2) - S_(m+1)^2
      = sum over i < j of u_i^m u_j^m (u_i-u_j)^2 > 0,
    u_j = 1/(x+jk).

The (0,1) term is positive. Finite termwise limits give x^r S_r(x) -> 1 at zero and -> p+1 at infinity, so both endpoint limits of R equal m/(m+1). Since every interior value is strictly smaller, R is neither nonincreasing nor nondecreasing. Adding n-2 common factors of derivative order m embeds this family into Conjecture 4.1 for every n >= 3.

The checker verifies the gap and strict deficit on 1,584 rational parameter choices, checks 168 multifactor embeddings and verifies both limits in nine rational-function cases. These finite checks are regression evidence; they do not certify the unbounded family or constitute a Lean proof. The general argument is the intended mathematical content of the subsequent formalization.

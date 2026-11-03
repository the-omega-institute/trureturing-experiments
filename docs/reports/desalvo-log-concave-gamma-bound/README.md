# DeSalvo Open Problem 8.2: Gamma obstruction

`check.py` uses only exact integer and `fractions.Fraction` arithmetic. It checks the proposed family

\[
p_0=\frac{M-1}{M^2+M-1},\qquad
p_{n+1}=\frac{M(M-1)M^{-n}}{M^2+M-1}\quad(M\ge 2),
\]

including normalization, positivity, log-concavity, the induced death rates, the bounded Stein solution for
\(A=\{0\}\), and the exact moments. It proves the exact lower bound

\[
\left(\frac{\Gamma}{\sigma}\right)^2\ge
\frac{M(M-1)^2}{2M^2-1}\ge \frac M8,
\]

so no finite universal constant can satisfy \(\Gamma\le C\sigma/w\) (take \(w=1\)). As a fidelity check,
the script also verifies the source's Bessel ratio identity exactly for \(k=3,\nu=5\). The output is stored
in `check.out`; the final line is `ALL_CHECKS_OK`.

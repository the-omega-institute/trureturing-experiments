# Hanna A281186 real-exponent positivity

Paul D. Hanna's sole COMMENTS conjecture in [OEIS A281186](https://oeis.org/A281186/internal), revision 10, asserts that the normalized formal solution of

\[
G=\exp\left(G^{-1}\int_0^xG^t\,dx\right)
\]

has entirely nonnegative coefficients exactly when the fixed real parameter satisfies \(t\ge 8/3\). The integral has constant coefficient zero, \(G(0)=1\), and \(G^t=\exp(t\log G)\).

Run `python3 check.py`. Python 3's standard library suffices. All evidence uses integers and `fractions.Fraction`; the output is in [check-output.txt](check-output.txt). Success ends with `ALL_CHECKS_OK` and exit code 0.

The checker reproduces all 20 published integer-normalized terms by literal formal fixed-point iteration and checks the source differential equation. It then evaluates an integer polynomial recurrence at degree 719 and independently recomputes the same integer by a Lagrange-inversion recurrence. The common negative integer gives

\[
[x^{720}]G_{8/3}=Q_{719}(0)/(3^{719}720!)<0.
\]

The SHA-256 of the signed integer's ASCII decimal representation, without a newline, is `effab943471a3e14c2c2c740a3375b3c5bf4195151ce57e51529f4d8f70f8acd`. The checker also verifies the rational exponential majorant used to compare the two real inverse obstructions uniformly for \(8/3\le t\le267/100\). Infinitely many negative coefficients throughout that interval requires the analytic inverse argument and Pringsheim's theorem; these arithmetic checks alone establish neither that analytic theorem nor a Lean settlement. No corrected sharp threshold is claimed.

The published terms and source equations are attributed to Paul D. Hanna and the OEIS Foundation, under the [OEIS end-user license, CC BY-SA 4.0](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement). The check program and its generated output are shared under CC BY-SA 4.0. No source article, command log or research transcript is included.

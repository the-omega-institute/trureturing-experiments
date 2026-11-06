# Hanna inverse-series coefficient checks

Run `python3 check.py` with Python 3; only the standard library is required. The
stored output ends in `ALL_CHECKS_OK`. All calculations use integers and rational
numbers. No floating-point value is evidence in this report.

The sources are Paul D. Hanna's [OEIS A258900](https://oeis.org/A258900/internal),
revision 18, and [OEIS A258924](https://oeis.org/A258924/internal), revision 11.
Their inverse-series definitions and displayed coefficient lists supply the
fidelity checks. V. Kotesovec's unnumbered FORMULA conjectures concern the exact
gamma constants in their coefficient asymptotics.

Two independent constructions are compared: Lagrange inversion of the specified
primitive, and the polynomial differential equation T'=1+T^p with T(0)=0,
followed by C=(1+T^p)^(-1/p). For p=4 and p=5 the program compares coefficients
through n=12, reproduces the 10 and 9 published sequence terms, respectively,
and reproduces five published S/C terms each. It checks F(S)=x, S'=C,
C^p+S^p=1 and T=S/C exactly. The same identities are checked through n=8 for
p=3,6,7,8.

The final checks verify the rational exponents in both proposed gamma
normalizations, conditional on the explicitly named gamma reflection and
positive elementary sine identities. They check the powers 17/6, 13/4 and
2^(14n) exactly. The program does not prove the reflection identities, analytic
continuation, a coefficient-transfer theorem or an infinite asymptotic limit.
The positive settlement requires a general analytic argument; finite checks
serve only to detect transcription and normalization errors.

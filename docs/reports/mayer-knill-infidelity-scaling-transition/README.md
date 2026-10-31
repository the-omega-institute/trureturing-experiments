# Mayer–Knill infidelity-scaling transition

`check.py` uses SymPy exact arithmetic to verify two pieces of the proposed
settlement.  First, for (d=2), ε = 1/10, and the depolarizing channel with
(p=2ε), the three probes (|0\rangle), (|1\rangle), and
((|0\rangle+|1\rangle)/\sqrt 2) each have input error ε and the process
fidelity is (1-3ε/2=17/20).  This reproduces the exact fidelity formula used
in the lower linear witness.

Second, the connected spanning family (|0\rangle) and
(v=(3/5,4/5)^T) admits the exact traceless non-scalar matrix
[
A=\begin{pmatrix}1&-3/2\\0&-1\end{pmatrix},
]
with (A|0\rangle=|0\rangle) and (Av=-v).  Thus the tested lines are
preserved while the nonorthogonality graph is connected, supplying the
square-root obstruction used by the proposed dichotomy.  The script uses
integers, rationals, and exact symbolic square roots only; its final line is
`ALL_CHECKS_OK`.

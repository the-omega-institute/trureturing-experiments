# Central-distance minimum size

`check.py` uses only Python exact integers and `fractions.Fraction`. It reproduces the source paper's bounded-band realization for a bicentral graph, verifies the round-robin factorization for every even order used by the proposed construction, and checks the proposed values

\[
m(1)=1,\qquad m(2)=2,\qquad m(k)=k\ (k\ge4\text{ even}),\qquad m(k)=k+1\ (k\ge3\text{ odd})
\]

for every `k` from 1 through 24. The output is stored in `output.txt`; the final line is `ALL_CHECKS_OK`.

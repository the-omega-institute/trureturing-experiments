# Higher logarithmic derivatives of chromatic polynomials

This exact check concerns Conjecture 2 of Fengming Dong, Jun Ge, Helin Gong,
Bo Ning, Zhangdong Ouyang and Eng Guan Tay, *Proving a conjecture on chromatic
polynomials by counting the number of acyclic orientations*, arXiv:1803.08658v3,
Section 6; Journal of Graph Theory 96 (2021), 343–360, DOI 10.1002/jgt.22617.
The conjecture asserts that every derivative of order at least two of
log((-1)^n P(G,x)) is strictly negative for x < 0.

Run with Python 3 and SymPy:

```sh
python3 check.py
```

`output.txt` is the program's standard output. It ends in `ALL_CHECKS_OK` on
success. The calculation uses integer and rational arithmetic exclusively.

The graph polynomial is independently computed from an edge list by summing
Whitney's spanning-subgraph expansion through frontier connectivity states.
The check agrees with direct coloring enumeration on all 75 graphs on one
through four labeled vertices, with one through four colors. As a result
stated in the source, it reproduces epsilon(K_n) = n - H_n, Section 1,
equation (7), for n = 1 through 6, using Proposition 1's logarithmic derivative.

For the theta graph consisting of three internally disjoint paths of edge
length a, the edge-list calculation agrees with the endpoint-coloring formula
for a = 2 through 9. For a = 7 the graph has 20 vertices and 21 edges. Its
fifth logarithmic derivative at -1 is exactly

```text
558835859111752710016430691303/5061397078105038614084902487500 > 0.
```

Three exact derivative computations agree: the fifth-derivative formula,
the rational Taylor recurrence obtained from Q' = Q (log Q)', and symbolic
differentiation of Q'/Q four times, where Q(t) = P(G,-t) for this even-order
graph. The change x = -t supplies the negative sign for the fifth derivative.

Identifying a single common vertex in m otherwise disjoint copies gives a
connected simple graph with 19m+1 vertices and 21m edges. Its chromatic
polynomial satisfies P_m(q) q^(m-1) = P_1(q)^m. The program independently
reconstructs these edge lists and their polynomials for m = 1 through 6,
checks connectivity, and verifies the fifth derivative equals m c + 24(m-1),
where c is the positive rational above. The algebraic decomposition
c + (m-1)(c+24) proves positivity for every integer m >= 1 once the standard
vertex-identification identity is established. Finite computation alone is
not a kernel proof of that identity or of the universal family statement.

The report supplies reproducible exact evidence for a preregistration draft.
It contains no Lean proof of the conjecture or its negation. The graph-counting
bridge, polynomial existence and uniqueness, and analytic derivative identities
remain obligations for subsequent formalization. Brown–Hickman's separate
root-location theorem is not an input to this program.

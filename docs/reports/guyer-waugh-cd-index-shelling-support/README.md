# Simplex shelling-component cd-index support

Exact checks for Question 3.16 of Dan Guyer and Alex Waugh, *Interpretations of the cd-index*, arXiv:2610.25491v1 (18 October 2026).

Run `python3 check.py`. The program uses only Python integer arithmetic and needs no external packages. Its expected output is in [output.txt](output.txt); a successful run ends with `ALL_CHECKS_OK`.

The first constructor counts flag chains in the finite face poset of the semisuspension of the first k facets of the n-simplex. Facets omit vertices 0 through k−1. It computes the flag h-vector by the Boolean Möbius transform, then exactly inverts c=a+b and d=ab+ba, rejecting any nonzero residual. Successive differences recover the individual shelling contributions. For all 28 initial segments with 1≤k≤n≤7 these contributions equal those computed by the separate implementation of the published Corollary 3.10 recurrence.

The recurrence constructor checks all degree-n cd words and all component indices for 1≤n≤14. It reproduces Corollaries 3.17 and 3.18 and checks the proposed full support classification: the initial component supports precisely words ending in c; the terminal component is zero; an interior component i supports a word u d c^t precisely when i+t+1≤n. It also checks the exact binomial coefficient expansion for every applicable trailing-c word through degree 14.

These are finite checks, not a proof for arbitrary dimension. The proposed general argument must establish the representation bridge, recurrence, nonnegativity, and terminal-tail induction. No claimed support pattern enters either polynomial constructor.

Source: https://arxiv.org/abs/2610.25491v1

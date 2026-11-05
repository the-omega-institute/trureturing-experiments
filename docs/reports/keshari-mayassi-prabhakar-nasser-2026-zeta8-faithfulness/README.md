# Multi-virtual twin group: zeta-eight faithfulness

Source: Vaibhav Keshari, Taher I. Mayassi, Madeti Prabhakar and Mohamad N. Nasser,
*Presentations and Representations of the Multi-Virtual Twin Group and Associated
Subgroups*, [arXiv:2605.13090v1](https://arxiv.org/abs/2605.13090v1), Question 4.3.

Run with Python 3 and SymPy (output generated with SymPy 1.14.0):

```sh
python3 docs/reports/keshari-mayassi-prabhakar-nasser-2026-zeta8-faithfulness/check.py
```

`check.out` is the exact output. The program uses rational functions over
Q(a,b,y0,y1,y2); b and every y parameter are nonzero. No floating point is used.
It reproduces the three matrix identities in Theorem 4.2(b) at a = -1, 0, 1,
and checks every instance of relations (11)-(19) for (k,n) = (2,3), (3,3), (3,4).
The last case also exercises all far-commutation families.

For k=3, n=3, put p = rho_1^0 rho_1^1 and q = rho_1^0 rho_1^2. Their images are
respectively diag(y1/y0,y0/y1,1) and diag(y2/y0,y0/y2,1), so [p,q] has identity
image for every a, b and every allowed y parameter. Mapping every s generator
to 1 and rho_i^alpha to the alpha-th free involution respects all 20 relation
instances. The image of [p,q] has reduced form (0,1,0,2,1,0,2,0), of length 8.
An independent action sending those involutions to (01), (12), (23) in S4 also
respects all relations and sends the word to (02)(13), certifying nonidentity.

The original k=2 specialization y1=-y0 is checked as well, with word
(rho_1^0 rho_1^1)^2 and reduced form (0,1,0,1). Its unrestricted matrix image is
printed to make its parameter restriction explicit. The k=3 commutator has no
such restriction and refutes the universal affirmative reading of Question 4.3.
It does not classify faithfulness for k=1 or k=2 with generic parameters, or for
other strand counts. These exact checks are preregistration evidence, not a Lean
proof or a formal settlement.

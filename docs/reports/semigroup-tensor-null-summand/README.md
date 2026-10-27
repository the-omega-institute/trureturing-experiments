# Tensor-power null summands for finite semigroups

Source: David He and Daniel Tubbenhauer, *Tensor powers of representations of
(diagram) monoids*, [arXiv:2508.04054v1](https://arxiv.org/abs/2508.04054v1),
Conjecture 4.2 and Example 4.3.

Run `python3 check.py`. The program uses only Python 3's standard library and
integer arithmetic. Its deterministic output is in `output.txt`.

The fidelity check reads the four-vertex fusion matrix printed in Example 4.3,
with vertex dimensions `(1,6,1,8)` and the third vertex the null module. It
reproduces 28 copies of the null module in the second tensor power and checks
the weighted dimensions and the formula `b(n)=6^n-7` for `2<=n<=12`. It checks
the consequences of the published matrix; it does not reconstruct the
partition-semigroup representation underlying that matrix.

The independent checks enumerate every nonempty family of singular binary
2 by 2 matrices (1023 families). They compare the proposed null multiplicity
formula with the pairing of the common kernel and the annihilator of the
common image, construct an invariant complement, and check the tensor bound
for exponents 1, 2 and 3. At exponent 1 they also exhaust every binary 2 by 2
idempotent with null image and null coimage, obtaining the maximum rank of a
null summand independently. Among these families, 92 are closed under
multiplication and 79 have no identity. An exact integer threshold is checked
on the full family at exponent 5. Binary field calculations remain valid
after extension to its algebraic closure.

The nilpotent Jordan checks enumerate coordinate chains of the tensor action
for `2<=d<=6` and `1<=n<=5`. The chains give actual invariant Jordan-block
decompositions, whose length-one block counts and total block counts match
`d^n-2(d-1)^n+(d-2)^n` and `d^n-(d-1)^n`, respectively. These coordinate
decompositions are valid over every field.

All checks pass with final line `ALL_CHECKS_OK`. These are bounded exact
computations supporting a preregistration; they do not prove the universal
conjecture or certify its literature status. The code is independently
authored; the fusion matrix data are attributed to the source above.

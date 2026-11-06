# Three-group character zeros

Run `python3 check.py` with Python 3. The program uses integers only. It represents
roots of unity in `Z[X]/(1+X^(q/3)+X^(2q/3))`, where `q` is a power of three.
`output.txt` is the complete output; the final line is `ALL_CHECKS_OK`.

The source is Alexander Moretó and Gabriel Navarro, *p-groups and zeros of
characters*, arXiv:2301.03917v2, Conjecture 2.6 and Theorem 2.7;
DOI [10.1007/s00013-023-01902-w](https://doi.org/10.1007/s00013-023-01902-w).

For every abelian group `U` of order 9 or 27, the program enumerates every
endomorphism in generator coordinates, retains automorphisms `alpha` with
`alpha^3=1`, and retains actions whose cyclic extensions have maximal class.
For every center-faithful linear character of each retained action it checks
that its induced degree-three trace vanishes inside `U` exactly on
`ker((alpha-1)^2) - ker(alpha-1)`. The inner product is computed exactly and
is 1. Every possible fixed twisting parameter `g^3=t` is included; a twist is
nonsplit relative to `U` exactly when `t` is outside the norm image.

For each of these concrete extensions, conjugacy classes are enumerated by
conjugation with generators. This reproduces Theorem 2.7 at `p=3`: the equality
`3^n-3^(n-1)+6` for zero elements corresponds to eight zero conjugacy classes.
The enumeration comprises 646 labelled actions, 1,938 extensions and 28
nonsplit twists. These counts are not counts of isomorphism classes, and do not
include groups without an abelian subgroup of index three.

Additional actions on the additive groups of `Z[omega]/(omega-1)^m`, for
`m=4,5,6`, test every center-faithful linear character. These cover selected
groups of orders 243, 729 and 2,187, with the same six internal zeros. Their
three twisting parameters have identical induced traces. These finite checks
do not establish the conjecture for arbitrary orders or reproduce O'Brien's
reported verification through order `3^10`.

A separate boundary check records why an unqualified literal interpretation
needs care: `C9` has class `1=2-1` and an abelian maximal subgroup of order three,
but the prescribed zero count is 12, exceeding the group order. Both arXiv v2
and the journal version omit an explicit nonabelian restriction in Conjecture
2.6 and do not supply a definition of maximal class. The boundary calculation
uses the convention “class n-1” without an additional nonabelian requirement;
it says nothing against the substantive conjecture on nonabelian groups.

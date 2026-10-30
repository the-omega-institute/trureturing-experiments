# Evolution-algebra graph balance under an annihilator shear

`python3 check.py` uses Python 3's standard library and exact integer/rational
arithmetic. Its stored output is `output.txt`; success ends in `ALL_CHECKS_OK`.

The fidelity check reproduces the single-arrow example in Section 1.1 of
Kaygorodov–Martín González–Viruel, *Evolution algebras: a survey of results and
a list of open problems*, [arXiv:2610.22248v1](https://arxiv.org/abs/2610.22248v1):
the algebra with `e1^2=e2`, `e2^2=0` has nilpotency index 3.

For the final question in Question 27, the program checks a family of dimension
`m+2`: `c_i^2=c_(i+1)`, `c_(m-1)^2=c_0+u`, `u^2=z`, `z^2=0`, with every mixed
product zero. Replacing `c_0` by `c_0+z` adds the arrow `c_(m-1) -> z`.
All 95 cases (`2 <= m <= 20`, over Q and the prime fields F2, F3, F5, F7)
have determinant-one basis transport, annihilator spanned by `z`, connected
graphs in both bases, and graph balances `m` and `1`. Every transported pair
product is computed. Matrix rows record squares; basis-change columns record
the new basis in old coordinates.

The graph check computes the gcd of fundamental closed-walk balances using
integer potentials on a spanning forest. This computes the all-length gcd:
tree potentials telescope on any closed walk, and each edge defect is itself
a fundamental closed-walk balance. Separate bounded enumeration cross-checks
all closed-walk states through length 7 for `m <= 6`, retaining the identity
of antiparallel arrows. The old graphs also satisfy the height-increment
condition modulo `m`; the new graphs contain a closed triangle of balance 1.

Path and balance definitions come from Alshatnawi–Costoya–Viruel,
[arXiv:2607.15154v1](https://arxiv.org/abs/2607.15154v1), Definitions 2.1, 2.2
and 2.4. The displayed formula counts a self-loop in both sets, giving zero;
the script preserves that literal formula. All counterexample graphs are
loop-free, so this ambiguity does not affect them.

These are finite exact checks and a source-fidelity check. They are not a
Lean proof of the universal question's negation or of the all-field family.

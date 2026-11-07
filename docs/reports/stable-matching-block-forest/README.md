# Stable matching block-forest realization

Source: Yoshiteru Ishida, *Prime factorisation of stable-matching instances:
uniqueness, simultaneous products, and an exact census*,
[arXiv:2610.05617v1](https://arxiv.org/abs/2610.05617v1), Section 7, Problem D.
The manuscript is dated 28 September 2026; arXiv v1 was submitted 4 October 2026.

Run with Python 3 and its standard library:

```sh
python3 check.py
```

`check.out` contains the exact output. Success ends in `ALL_CHECKS_OK` with exit 0.
All evidence uses integer and finite-set arithmetic.

The fidelity check enumerates every balanced strict complete preference profile
at sizes one, two and three. Two independent block detectors, using top-set
prefixes and literal rank inequalities respectively, agree throughout. Exact
covers of both sides reproduce Theorem 9 and Table 3: 14 prime profiles out of
16 at size two, and 45,768 prime profiles out of 46,656 at size three; the other
size-three finest partition counts are 504 for 2+1 and 384 for 1+1+1.

The construction check generates every unordered positive-size rooted tree with
root size n, smaller children and child-size sum at most the parent size, for
1 ≤ n ≤ 7. The shape counts agree with Problem D: 1, 3, 10, 40, 170, 785, 3,770.
Each prescribed child receives its recursive preference lists. Uncovered agents
are separate singleton cells; cyclic tails precede the uncovered proposer's own
receiver and follow every other cell's internal lists. Receivers rank their
own indexed proposer first. Four variants reverse the cell cycle order and/or
the enumeration inside cells. Every run checks strict completeness and equality
of the entire block family with the prescribed family, including the absence of
blocks with different proposer and receiver index sets. There are 19,116 runs;
the variants can coincide. The prime-compatible counts also equal 2 f(n−1) for
2 ≤ n ≤ 7. The supplied size-five nested and side-by-side example is reproduced.

This independently written checker imports no upstream code. Its finite results
support the proposed general construction; they do not prove universal
realizability or provide a Lean proof.

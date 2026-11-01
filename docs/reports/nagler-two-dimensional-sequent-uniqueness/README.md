# Nagler two-dimensional sequent uniqueness

`check.py` is an exact finite check for the preregistered refutation of
Conjecture 1 in Sophie Nagler, *Inference-Behaviour Semantics for All*
Connectives in Two-Dimensional Sequent Calculi*, arXiv:2607.19419v1.

It reproduces the source's stated 376-member `PC` count, records that the
parenthesized PU bullet counts sum to 152 although the source prints 150, and
verifies that the candidate `⟨α, δcd⁻ 1B, 1, N⁰:N⁰⟩` is in the PC family and outside the listed PU
families. It constructs the two three-rule bridges by substituting `c(A)` into the
independent `B` guard, and checks both directions of all four contextual Cut
replacements for finite multisets of size at most two.  All arithmetic is
integer/tuple arithmetic; there is no floating-point evidence.

Run:

```sh
python3 check.py > check.out
```


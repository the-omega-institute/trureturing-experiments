# Cassaigne Open problem 2

`check.py` uses Python integers, `fractions.Fraction`, and SymPy's exact
`sqrt(5)`. It reproduces the paper's occurrence of `(aba)^3` at position 5 in
the Fibonacci word and checks the proposed family `v_k` obtained by replacing
each Fibonacci letter by a block of `k` equal letters. For `k = 2, ..., 6`, the
long generated prefix contains all four bigrams and satisfies
`p_{v_k}(n) = n + 2k - 1` for every tested `n >= 2k + 1`.

These are finite checks supporting the preregistration. They do not claim to
prove the asymptotic limits or the infinite-family theorem.

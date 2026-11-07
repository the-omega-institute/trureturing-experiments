# Elek magic-word converse: exact preregistration check

`check.py` uses only exact integer arithmetic. It reproduces the paper's
S4 example (the word `(s1,s3,s2,s1,s3)` has one positional subword with
product `s2`), then checks the proposed B2 witness `v=s`, `w=sts`.
It enumerates all reduced words and positional subwords, verifies the four
interval elements form a Boolean lattice, and counts the two subwords for
`s`. The same dihedral normal-form calculation is checked for every
`I_2(m)`, `m=4,...,100`; the accompanying preregistration gives the symbolic
all-`m >= 4` argument. No floating-point arithmetic is used.

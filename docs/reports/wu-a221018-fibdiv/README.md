# Wu’s OEIS A221018 decimal-residue conjecture

Run `python3 check.py` with Python 3; the standard library is sufficient.
`output.txt` is the exact output. All evidence uses integer arithmetic.

Source: [OEIS A221018](https://oeis.org/A221018), revision 35, first FORMULA
conjecture attributed to Chai Wah Wu. The included `b221018.txt` is Wu’s table
of the first 10000 terms, obtained from https://oeis.org/A221018/b221018.txt.
The OEIS data are attributed to the OEIS Foundation and the named contributors
under [CC BY-SA 4.0](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement).

The program reproduces every member and gap for 2 ≤ k ≤ 500000, including the
17 terms in the source Python range, and the source’s count of 81 non-8 final
digits among the first 7777 terms. It checks all 10000 listed memberships and
the conjectured final digits. A direct integer recurrence cross-checks every
index through 20000. The small-period rank ingredients are checked explicitly.

The final-digit histogram of the 10000 listed terms is {0: 86, 2: 1, 8: 9913}.
The program also checks k=14930352 (term 126) by fast doubling, independent
matrix powering and 14930352 direct modular recurrence steps. This term is
12 modulo 30, contradicting the stronger comment in OEIS A225219, while it
is 2 modulo 10 and satisfies A221018’s conjecture.

These are finite checks and source-fidelity evidence. They do not establish
the unbounded conjecture. No Lean proof is included.

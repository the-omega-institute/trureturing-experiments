# Six-letter three-ball intersection checks

The source is Junnila, Laihonen, Lehtilä and Padavu Devaraj, *The Size of the Intersection of q-ary Hamming Balls*, arXiv:2606.09158v2, Section III-A2 and Section IV.

Run `python3 check.py` with Python 3 and SymPy. All evidence uses integers or exact symbolic rational arithmetic. The output is stored in `check-output.txt`; success ends with `ALL_CHECKS_OK`.

The program verifies the polynomial factorization, derives the universal next-radius coefficient symbolically through a fixed Laurent multinomial sum, checks integer sparse enumerators for k=1 through 15 and three common-coordinate lengths, and enumerates actual words for the smallest odd-distance configurations. For alphabets of sizes 3, 4, 5 and 7 it reproduces the source’s central-word counts q and 6. For alphabet size 6 it obtains 6 versus 6 at radius 2 and 189 versus 246 at radius 3 in dimension 4.

These checks support the preregistered paper-level argument. They do not prove the global extremal reduction or establish Lean kernel verification. The proposed positive settlement uses an argument uniform in distance and radius, with these computations as supporting checks.

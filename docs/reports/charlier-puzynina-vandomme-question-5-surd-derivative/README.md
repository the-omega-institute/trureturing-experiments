# SURD derivative: exact finite controls

The source is É. Charlier, S. Puzynina and É. Vandomme, *Recurrence along directions in multidimensional words*, arXiv:1907.00192v2, Section 8, Question 5. The proposed counterexamples are the binary words

    w_m(i,j) = 1 iff (i+j) mod m = m-1, for m >= 4.

Run `python3 check.py`; its complete output is in [check.out](check.out). The program uses only exact Python integers and the standard library. It independently reproduces the source's four Thue–Morse return words to 011 and the derivative prefix 12341243123431241234124, using the source's one-based labels in that example. The multidimensional derivative uses labels starting at zero, separately on each primitive ray.

The remaining controls check 242,448 block-periodicity instances for m=4,...,12; 1,310 primitive rays in the three prescribed coordinate-sum residue classes; and 63 arithmetic constructions with m=4,...,12 and N=1,...,7, yielding 252 exact-gcd and adjacent derivative-block checks. The composite-modulus control uses 551=19*29. The largest constructed coordinate has 2,460 bits.

A derivative value at a large lattice point is computed by gcd-normalizing that point and enumerating consecutive returns in two full periods of its sampled input. This checks that the return-word cycle repeats and then indexes it by the radial multiplier. It does not enumerate an exponentially large prefix or assume the three proposed ray formulas.

Every tested obstruction has initial shape-(2,1) derivative block 00 and N consecutive sampled blocks 01. These finite controls do not prove non-SURD, whose quantifier requires every N, and do not certify any Lean theorem. The symbolic family argument and literal binding definitions belong in the preregistration draft. The check source and report are project-authored; no third-party program or paper text is vendored.

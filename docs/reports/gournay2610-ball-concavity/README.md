# Even-ball growth in Boolean free products

For Antoine Gournay, *Entropy, isoperimetry and the total variation profile of random walks*, arXiv:2610.16834v1, Question 6.9.

Run `python3 check.py` with Python 3.12 or newer. Only the Python standard library is required. `output.txt` contains the exact output. Every arithmetic operation uses integers or `fractions.Fraction`; no floating point is used.

The source-fidelity check reproduces Section 5: the support of a lazy walk at time n is the closed word ball, and the even return probability is the squared L2 norm. The implementation uses probability 1/2 at the identity and 1/(4k) at each coordinate involution of (C2)^k * (C2)^k.

Independent breadth-first enumeration of reduced bitmask syllables checks the weighted normal-form recurrence for k=1,2,3. Enumeration of positive compositions checks the recurrence at radii 4,6,8,10,12,14 for k=2,7,8,12,13,14,20. Exact polynomial arithmetic reconstructs b4,b6,b8, the identity b4*b8-b6^2=k*Q(k)/1209600, and every positive coefficient of Q(13+t). This certificate gives a violation for every integer k>=13 at central radius 6. The k=8 violation at central radius 12 is checked separately. The classical C2*C3 even-ball formula supplies a control.

These are finite exact arithmetic and symbolic identity checks, conditional on the mathematical identification of normal-form weights with word length. They are not Lean proofs of the group metric, exponential growth, or the universal family statement. The infinitely-many-radii pole argument is not certified by this program.

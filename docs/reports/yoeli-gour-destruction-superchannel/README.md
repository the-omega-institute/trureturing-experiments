# Destruction-covariant channels and fixed superchannels

Source: Goni Yoeli and Gilad Gour, *Instability as a Quantum Resource*,
[arXiv:2602.18323v1](https://arxiv.org/abs/2602.18323v1), Open directions,
question after Eq. (29).

Run `python3 check.py` with SymPy installed. All computations use exact
integers, rational numbers and Gaussian rational numbers. No numerical
tolerances or floating-point evidence are used. The stored output is
`check-output.txt`; a successful run ends with `ALL_CHECKS_OK` and exits 0.

The source fidelity check verifies the complete dephaser's Kraus formula,
idempotence and full-rank fixed state in dimensions 1, 2 and 3. It reproduces
Eq. (29) on every matrix-unit linear map: pre- and post-dephasing induces
complete dephasing of the Choi matrix. Its fixed channels are classical;
the identity channel is destruction-covariant but not fixed when d >= 2.

For d = 2, 3 and 4, the program constructs all permutation unitaries and
their single-coordinate phase-i variants, checks their DIO property, and
produces every vectorized matrix unit by differences. Their exact nonzero
overlap graphs have respectively 2, 1 and 1 components. A separate exact
linear-algebra computation gives joint commutant dimensions 2 and 1 at
d = 2 and d = 3. The overlap identities used by the general argument are
checked throughout these families. A rational unitary with a 3/5, 4/5
rotation block is outside DIO in each tested dimension.

The qubit check verifies, on the entire space of linear maps, that the
physical involution twirl N -> (N + C N C)/2 is idempotent and fixes exactly
DIO. Here C is conjugation by diag(1,-1). The attempted channel involution
2 Delta - id is not completely positive at d = 3: its Choi quadratic form
on vec(I) equals -3.

These are finite checks supporting a dimension-uniform ordinary argument.
They do not certify the infinite family or constitute a Lean proof. The
proposed obstruction uses the fact that a unital CP map fixing a projection
has Kraus operators commuting with that projection, followed by the
spanning and connected-overlap argument. Only the eventual formal proof
can certify that argument in the repository's formal system.

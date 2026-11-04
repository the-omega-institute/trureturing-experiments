# Positive pure branch on [11/5, 9/4]

This exact certificate supports FourPath Continuation II, Section 50. It uses
the original ordinary transpose, one whole feedback scalar and the original
joint distance. The mathematical suppliers are the original FourPath volume,
Continuation I Section 37 and Continuation II Sections 42, 43, 46 and 47 at
source pin `5fd6a884f9173f778cfa040217d7d51f748ac9c6`. Those three volumes are
byte-identical at the theory delivery baseline
`a6e2c138042b4ba25416edf5693b51d19678ed3e`.

For every level xi in the closed interval [11/5, 9/4], every original positive
normalized rank-one effect and every original U in U(4) whose positive distance
branch is active at true joint distance xi, the result is g(E,U) >= V(xi).
Here V is only the selected positive rank-one geometric value of Theorem 37.5.
The selected determinant is positive. If the actual population slice has
negative positive-row determinant, then

    g(E,U) - V(xi) > 1196031/1600000000 > 7/10000.

The certificate at (xi,q)=(9/4,1/4) has omega < -19/1000, so the excluded
negative-determinant slices are nonempty. The Section 47 gain estimate and
Section 42 original-unitary transport are reused; this program certifies the
continuous scalar inequalities consumed by the new interval comparison.

Run with Python 3.12 or newer, from any directory:

```sh
python3 -B /path/to/certify_interval.py
```

The program deterministically replaces its sibling `certificate.json`.
All proof arithmetic uses standard-library `Fraction` and integer `isqrt`;
there is no floating-point sign test. Let S=2^96. At each rational root endpoint
r >= 0, n=isqrt(floor(r*S^2)); the lower root is n/S and the upper root is
(n+1)/S, except that an exact square uses n/S for both. Endpoint-square
containment is asserted. Addition, subtraction, four-corner multiplication,
zero-excluding reciprocal, nonnegative squaring and square root are natural
interval extensions. The constants are h=sqrt(1/2), c=sqrt((1+h)/2),
s=sqrt((1-h)/2), k=c-h, with exact identities h^2=1/2 and s^2=(1-h)/2.
The formulas for omega and the fixed-xi derivative are those in Section 50.

| Certificate | Function and closed domain | Required bound | Pieces | Outward aggregate enclosure |
| --- | --- | --- | --- | --- |
| C1 | omega(9/4,q), q in [0,1/10] | >= 1/10000 | 10 | [313809/250000000, 23008109/200000000] |
| C2 | omega(9/4,q), q in [3/8,1] | >= 1/10000 | 27 | [129303/1000000000, 19833451/31250000] |
| C3 | omega(9/4,q), q in [1/10,3/8] | >= -21/1000 | 132 | [-20992343/1000000000, 4663773/250000000] |
| C4 | G'_xi(9/20), xi in [11/5,9/4] | <= -1/10 | 1 | [-453826117/1000000000, -24250863/125000000] |

Each domain is recursively bisected at rational midpoints until the interval
image of the entire closed leaf satisfies its threshold. The manifest stores
every rational leaf endpoint and its full rational image, not just these
aggregate bounds. Coverage assertions check the original two endpoints,
strictly positive leaf lengths, equality of adjacent endpoints and exact total
length. Section 50 also supplies a compact encoding of every partition.
Aggregates round the minimum lower endpoint down and the maximum upper endpoint
up to multiples of 10^-9. Exact manifest constants include the curvature,
positive denominators' domain, gain bound, strict margin and inactive-sign gap.

| Artifact | Bytes | SHA256 |
| --- | --- | --- |
| `certify_interval.py` | 6655 | `edb176c129bce7c889c0bd64d85eb87e899acf74503435fa6f429f280e210a95` |
| `certificate.json` | 186021 | `661b9582e4425abb56e4fc17e73d7aa68a2b5970233c0a95f1109c33a85f4046` |

This is exact rational certificate evidence for an ordinary mathematical proof,
with no Lean/kernel, freeze, coverage or final review claim. It does not settle
the negative-active branch, mixed-input pointwise comparisons, the original
entire range, a prefix from zero, first loss, alpha/beta, full equality/lift
classification or the continuing AURIC goal. Existing selected pairs supply
attainment. No literature priority is claimed.

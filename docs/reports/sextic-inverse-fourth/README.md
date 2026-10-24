# Actual Eisenstein inverse fourth moments

These programs evaluate one fixed component of the sextic Hecke family over
`Z[omega]`, where `omega^2+omega+1=0`. They enumerate every good ideal in each
specified norm annulus and all six unit multiples of each row. Prime ideals
above 2 and 3 are excluded. Rows are sixth-power-free; inverse columns are
squarefree. The fixed ray twist is trivial. This is an actual arithmetic family,
not independent random characters or a prime-row surrogate.

The sextic symbol `chi_d(u)` uses the residue of
`u^((Norm(p)-1)/6)` at each prime ideal `p | d`, with its multiplicity in `d`.
The root labels are the powers of `zeta_6=1+omega`. Split and inert rational
primes use their actual quotient fields. If a row meets a denominator prime,
the symbol is zero, including at every positive power of that character.

For each integer scale `U`, put `D=sqrt(U)` and evaluate

$$
M(u)=D^{-1/2}\sum_d\mu_K(d)V(Nd/D)\chi_d(u),\qquad
\mathcal T=\sum_u\Omega(Nu/U)|M(u)|^4.
$$

`V` is the smooth bump centered at 1 with radius 1/2, and `Omega` is centered
at 3/2 with radius 1/2. The normalized bump at `|t|<1` is
`exp(1-1/(1-t^2))`, and is zero outside. The programs also evaluate the same
second and fourth moments after replacing `mu_K` by `|mu_K|`, the ideal
Liouville function, or the constant 1. All four use the same rows and zeros;
the ordinary and Liouville comparisons include nonsquarefree columns.

The literal inverse square is grouped by its parity core:

$$
M(u)^2=\sum_a g_D(a)\chi_a(u),\qquad
H_2=\sum_s\sum_u\Omega(Nu/U)
 \left|\sum_{r_2(a)=s}g_D(a)\chi_a(u)\right|^2,
\qquad E_6=\mathcal T-H_2.
$$

`E6` is signed. `H2` includes its within-block interference, and is not
replaced by the diagonal. A second computation retains exactly the terms
`gcd(d1,d2)=1` in `M^2`. Its coefficient columns are squarefree, so its parity
blocks are single columns. The report records that head's energy, diagonal
and signed difference, as well as the tail energy and full head/tail
interference. It never infers the full moment from the head alone.

Six-unit averaging is computed by finite Fourier orthogonality. An independent
scalar checker instead visits every unit and every ordered factor pair. The
reported reconstruction error compares the grouped inverse square against
`M(u)^2` on every enumerated row. All eligible rows are enumerated; there is no
sampling of the arithmetic family.

## Run

Requires Python 3.9 or newer and NumPy. From the experiments repository root:

```sh
python3 -B docs/reports/sextic-inverse-fourth/test_eisenstein.py
python3 -B docs/reports/sextic-inverse-fourth/test_moment.py
python3 -B docs/reports/sextic-inverse-fourth/moment.py --scales 256 1024 4096 --out /tmp/sextic-inverse-fourth-pilot.json
python3 -B docs/reports/sextic-inverse-fourth/moment.py --scales 16384 65536 --out /tmp/sextic-inverse-fourth-larger.json
```

The output path is explicit and is overwritten if it already exists. Source
hashes, dependency versions, row and column counts and computation times are
included in each report. `--chunk` controls temporary memory only.

## Interpretation and sources

Character exponents, ideal factorization, unit actions and collision zeros are
computed with exact integers. Smooth profiles and the reported moments use
binary64 arithmetic, with independent scalar comparisons; the outputs are not
interval certificates. Finite ratios or slopes do not establish an asymptotic
power saving, a uniform inverse-fourth estimate, Robin's inequality or RH.

The mathematical object and remaining analytic target are specified in
[the analytic source note](https://github.com/the-omega-institute/trureturing/blob/6321bef543a32f29b449551645cbf2e391d2848d/Library/Analytic/liu2026quasirh.md),
particularly its actual inverse coefficients and order-six correlation. The
sextic-residue, primary-generator and unit conventions follow
[Collet's pinned foundational exposition](https://github.com/rohanarun/quasi-rh-hillclimb/blob/594fc4838e9c92d5a16bedf270d07d1f51161be7/submissions/hailey-collet-1252270103/paper/foundations.tex),
Sections 1.2 and 2.1. The arithmetic algorithms and experiment code here are
independent implementations; no external program source is vendored.

## Retained finite results

The complete outputs are [pilot.json](pilot.json) and [larger.json](larger.json).
Values below are rounded displays of their binary64 measurements.

| U | Physical rows | Inverse fourth moment | H2 | E6 | Shared-gcd-one E6 |
|---|---:|---:|---:|---:|---:|
| 256 | 480 | 13.619116 | 14.291846 | -0.672730 | -0.736723 |
| 1024 | 1872 | 38.246197 | 39.414064 | -1.167868 | -0.885732 |
| 4096 | 7386 | 117.323817 | 117.964805 | -0.640988 | -0.099621 |
| 16384 | 29760 | 655.537800 | 661.323346 | -5.785546 | -7.488157 |
| 65536 | 118824 | 2686.143572 | 2640.250156 | 45.893416 | 42.509043 |

The order-six cross-core term changes sign across these finite scales, including
in the shared-gcd-one head. Its observed negative sign at the smaller scales
cannot be used as a uniform sign rule. At U=65536 the Möbius fourth moment is
2686.143572, while replacing its signs by positive coefficients on the same
squarefree support gives 2649.062665. Thus the measured sign contribution is
not uniformly a reduction of the positive-coefficient moment either.

The largest relative inverse-square reconstruction error is below 1.4e-15.
The arithmetic checker covers 14 small prime-ideal fields, including inert
fields of sizes 25 and 121; the independent scalar checker makes 56 full-unit
moment comparisons. These checks support the implementation. The tabulated
signs are numerical observations, not interval-certified counterexamples or
asymptotic estimates. The missing analytic input still has to bound the
positive cross-core contribution with its actual Möbius coefficients.


## Product-squarepart correlation

`squarepart.py` resolves the full original order-six correlation by

$$
ab=p^2k,\qquad k=r_2(ab)\text{ squarefree},\qquad P=Np.
$$

It retains precisely the ordered unequal-parity pairs with $P^4<U$,
and also computes their complement. The strict threshold uses integers.
All original sixth-power-free rows, all six unit multiples, actual
Möbius convolution coefficients and collision zeros remain. The program
imports the existing ideal engine and smooth coefficient plan; it does
not replace the family by independent random characters.

The flat calculation assigns weight one to the selected inverse and row
annuli. Since the norm grids are finite, these values can be realized
by smooth profiles supported in the same open annuli, with transitions
between adjacent sampled norms. Such a profile is fixed for that finite
example; the flat runs at different scales are not asserted to come
from one common fixed smooth profile or to have uniform derivative costs.

Flat arithmetic uses the exact representation
$\zeta_6^j=A_j+B_j\zeta_6$, $\zeta_6^2=\zeta_6-1$.
The normalized fourth moments have denominator $U$, so every flat
quantity is a rational number obtained by integer arithmetic. An explicit
bound rejects requests that could overflow the integer array accumulators.
The second calculation uses exactly the fixed smooth bumps from `moment.py`;
these weights and sums use binary64 and are not interval enclosures.

The retained [squarepart.json](squarepart.json) evaluates the preselected
scales $64,128,256,512,1024,2048,4096$. In particular:

| U | Flat small-P correlation, exact | Fixed bump small-P correlation |
|---|---:|---:|
| 128 | $-9/4$ | -0.1603369703 |
| 256 | $-135/32$ | -1.1848357344 |
| 1024 | $-135/32$ | -1.2508545208 |
| 2048 | $747/64$ | -1.2448873346 |
| 4096 | $7041/512$ | 0.2490045477 |

The exact positive value at $U=2048$ prevents a blanket all-scales
nonpositive-sign assertion for this small-$P$ contraction over admissible
finite smooth profiles. It does not disprove eventual nonpositivity for
a specified fixed profile. The fixed-bump measurements also have both
signs, but those signs are numerical observations rather than interval
certificates. None of these outputs proves an asymptotic power bound,
a full inverse-fourth improvement, Robin's inequality or RH.

`test_squarepart.py` independently checks the positive flat example by
visiting the four original ordered inverse-factor slots and every
physical unit row, without the grouped inverse plan, pair plan or unit
orthogonality used in the main program. It also compares the new smooth
calculation's full moment and cubic block against the retained original
baseline and checks the exact signed decomposition
$\mathcal T=H_2+E_{6,\mathrm{small}\ P}+E_{6,\mathrm{large}\ P}$.

```sh
python3 -B docs/reports/sextic-inverse-fourth/test_squarepart.py
python3 -B docs/reports/sextic-inverse-fourth/squarepart.py --scales 64 128 256 512 1024 2048 4096 --out docs/reports/sextic-inverse-fourth/squarepart.json
```

The corresponding analytic target and the already controlled large-$P$
sector are [C16–C20 in the owning source](https://github.com/the-omega-institute/trureturing/blob/3abf0112f94615ae920928b3b01dcca20be9a113/Library/Analytic/liu2026quasirh.md#localizing-the-remaining-contraction-by-the-product-squarepart).
The remaining estimate must control positive contributions with the
actual coefficients; the Möbius sign cannot simply be assumed to pay it.

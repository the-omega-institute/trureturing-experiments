# Singleton-survivor probability checks

Source: K. P. O'Keeffe, *Arithmetic of the sync basin for pulse-coupled
oscillators*, arXiv:2609.01668v2, Supplementary Material §1, subsection
“The singleton-survivor volume for gamma <= 0”, equation `eq:sm:vsing`.
The preregistration target is N >= 2 and -1 <= gamma <= 0.

Run `python3 check.py` with Python 3, SymPy and NumPy. The complete output is
in [output.txt](output.txt); its last line is `ALL_CHECKS_OK`.

All computations use integers, rational arithmetic, or exact symbolic algebra.
The program reproduces the source's N=2 probability by independently enumerating
the original event polytope's vertices at gamma = -1/4, -1/2, -1. It also:

- Enumerates the original event polytopes through N=5 at those three parameters,
  checking their full vertex sets against the proposed weighted simplices.
- Checks inverse simplex vertices against every original event and band facet,
  and computes exact determinant volumes, through N=10 at gamma = -1/4, -1/2,
  -1, -4. The last parameter tests a range outside the preregistration target.
- Checks the gamma=0 staircase simplices through N=10 and the unit-triangular
  coordinate Jacobians symbolically in q.
- Compares the homogeneous event criterion with direct rational voltage
  simulation on 1,080 random configurations, and follows explicitly constructed
  interior survivors for 5N events.

The Monte Carlo uses independent uniform rational-grid initial voltages r/2^30,
where 1 <= r < 2^30, and integer homogeneous firing maps. Arithmetic overflow is
excluded by an explicit signed-int64 bound. The finite-firing criterion is the
source's characterization of eternal survival. These samples do not establish
the continuous probability identity. Their fixed discrepancy flag is six
standard deviations, expressed without square roots as a rational squared error.

| N | gamma | trials | survivors | observed frequency | formula |
|---|---|---:|---:|---|---|
| 6 | -1/2 | 10000000 | 570 | 57/1000000 | 27/495040 |
| 7 | -1/2 | 10000000 | 19 | 19/10000000 | 27/7235200 |
| 7 | -1/2 | 100000000 | 347 | 347/100000000 | 27/7235200 |

The two N=7 readings use the same seed; the 10-million sample is the prefix of
the 100-million sample. Their squared standardized errors are respectively
219556487168/24418708875 and 224168159296/122093544375. The N=6 squared error is
205801472/185629875. The lower N=7 pilot count is retained explicitly.

The check supplies finite exact consistency evidence and statistical readings.
It does not prove the universal statement, forward invariance for arbitrary N,
or any Lean theorem. No source text or third-party implementation is copied
into the program.

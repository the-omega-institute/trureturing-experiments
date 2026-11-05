# Shifted Stirling repair factors

Source: Piotr Miska and Tom Ward, *Stirling numbers and periodic points*,
arXiv:2102.07561v2; Acta Arithmetica 201 (2021), 421–435,
DOI 10.4064/aa210309-9-8. Target: Conjectures 1, 2, and only Conjecture 4(d).

Run `python3 check.py` with Python 3. The program uses only standard-library
integer arithmetic and `fractions.Fraction`. Its output is in `output.txt`.
Successful completion exits 0 and ends with `ALL_CHECKS_OK`.

The fidelity checks reproduce the source's second-kind sequence at k=3 and
all 40 hypothetical repair factors in Table 1 from denominator lcms at n≤3000.
They also verify nonnegative Möbius sums and the source factorial denominator
bound on the 120,000 pairs 1≤k≤40, 1≤n≤3000.

The target checks test the radical lower divisor and lcm upper divisor for
those truncated lcms, 19 sharp prime-power cases p^j+1≤40 and their generating
function residue formula through n=3000, and 612 polynomial quotient
congruences for 3≤k≤18, 1≤a,n≤3, and every prime p<k. The polynomial check
uses exact binary exponentiation modulo both p^a and the monic polynomial
∏(X−i), 0≤i≤k. It also checks the annihilation identity using an independently
computed Stirling recurrence.

A truncated denominator lcm is a lower divisor of the infinite repair factor.
These finite computations do not certify the infinite repair factor or prove
the universal conjectures. The source's additional tests through n=50000 are
not rerun here. The proposed settlement uses a general argument.

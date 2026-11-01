# Hardin A204614: width-seven ordered-minor recurrence

Run `python3 check.py > output.txt` in this directory. The program uses Python
integer arithmetic, exits zero on success, and prints `ALL_CHECKS_OK` last. It
also regenerates `certificate.json`. No network access or dependencies are needed.

The source is R. H. Hardin's [OEIS A204614](https://oeis.org/A204614), revision
9 (2025-10-03), originally submitted 2012-01-17. Its sole FORMULA asserts an
empirical order-43 recurrence for every n > 47. The exact entry, all 14 displayed
terms and all 210 linked b-file terms are in `source-data.json`. Source URLs and
the [OEIS license](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement)
are recorded there. Source data are attributed to Hardin and OEIS contributors
under the OEIS end-user terms, including CC BY-SA 4.0 for the main collection.

An array has n+1 binary rows of width seven. Adjacent 2 by 2 determinants are
integer determinants; their six entries in each minor row are nondecreasing,
and successive minor rows are coordinatewise nondecreasing. The program builds
the state (last row, last minor row) directly from all 128 squared row pairs.
It finds 265 states, 2,470 initial pairs and 3,051 directed edges. It verifies
every state, initial multiplicity and possible edge independently using integer
bit encodings, and directly counts 62,031 valid three-row arrays.

With source rows and target columns, let T be the 0/1 transition matrix,
alpha the row vector of initial multiplicities and 1 the all-ones column.
Deleting/appending the last row gives the faithful counting correspondence
a(n) = alpha T^(n-1) 1 for every n >= 1: the initial multiplicity records all
first-row choices, and every appended row is determined by its target state.
This explains why no array multiplicities are lost by the state compression.

Let P(X) = X^43 - sum(c_i X^(43-i)), using the source coefficients. The exact
program computes b = T^4 P(T) 1 and verifies alpha T^k b = 0 for all
k = 0,...,264. Independent prefix-vector counts agree with the suffix-vector
calculation, all 210 source terms, and all 265 recurrence residuals at
n = 48,...,312. At n = 47 the residual is -79321900236876546048.

The unbounded mathematical argument consumes these finite facts. The monic
characteristic polynomial Q of the 265 by 265 integer matrix T has degree 265.
Cayley-Hamilton gives Q(T)=0, hence an order-265 homogeneous recurrence for
r(k)=alpha T^k b. Its first 265 values vanish, so strong induction gives r(k)=0
for every k >= 0, exactly the published recurrence for every n > 47. There is
no need to compute the coefficients of Q. This general argument, together with
the array correspondence, supplies the universal step; finite testing alone
does not. The checker does not execute an infinite induction or a Lean proof.

`certificate.json` contains all states, edges, initial weights, b, 312 counts,
the 43 source coefficients and 265 zero residuals. It is generated from the
literal array conditions, not from a recurrence fit or the search seat's data.
The intended Lean implementation must prove the correspondence, check all
finite identities in the kernel, and reuse Cayley-Hamilton. It should use
balanced trees for certificate lookup and separate per-item `decide +kernel`
checks, rather than one large reducible conjunction. This directory establishes
exact computational evidence and a mathematical route, not kernel verification.

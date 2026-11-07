# Sibirmovsky A276436: exact support checks

Source: [OEIS A276436](https://oeis.org/A276436), revision 35, and
[Yuriy Sibirmovsky's Math.SE question 1911465](https://math.stackexchange.com/questions/1911465).
The conjecture says that the coefficientwise limit of the nested polynomials
has precisely the zero coefficients indexed by 0, 3, 5 and 8.

Run with Python 3 using only the standard library:

```sh
python3 check.py
```

`check.out` contains the result. All arithmetic used as evidence is integer
arithmetic. The script:

- reproduces every coefficient in the OEIS DATA field and the 51-term prefix
  displayed in the original question for the polynomial f_10;
- verifies the source's coefficient-stabilization bound for f_1 through f_16,
  within degree 100;
- computes the finite nested polynomials modulo X^101 by direct multiplication,
  obtaining zeros exactly at 0, 3, 5 and 8;
- independently counts internal-level profiles with binomial weights and checks
  all 101 coefficients against the polynomial calculation;
- constructs profiles for every nonexceptional index through 100,000 and checks
  the branching capacities, leaf count and edge count exactly.

At internal level d a profile t_d has (d+1)t_d available child slots.
It must satisfy t_1=1 and t_(d+1) <= (d+1)t_d. Its leaf count is
1 + sum(d*t_d). Write n = 1 + h(h+1)/2 + r, with 0 <= r <= h.
For r=0 use a chain. For r>=2 double t_r. For r=1 and h>=4 use a
chain through h-1 and double t_2 and t_(h-1). The exceptional small cases
are precisely the four rejected indices.

These are finite reproducibility checks for a proposed general proof.
They do not establish the unbounded conjecture and are not Lean verification.

The 50 OEIS coefficients embedded in the script are attributed to Yuriy
Sibirmovsky under the [OEIS CC BY-SA 4.0 terms](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement).
The Math.SE post is attributed to Yuriy Sibirmovsky under its CC BY-SA 3.0
content license. The algorithms in this directory are independently written.

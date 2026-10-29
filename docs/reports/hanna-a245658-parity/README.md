# Hanna A245658 parity: exact finite consistency checks

Sources: [OEIS A245658, revision 24](https://oeis.org/A245658) and
[OEIS A247031, revision 12](https://oeis.org/A247031), both 2014.
The check reproduces all 15 displayed coefficients of both entries from the
binomial-cube triangular transform and normalized cube-root recurrence.

Run with Python 3 (standard library only):

```sh
python3 check.py > check.out
```

All arithmetic uses Python integers. The program checks the defining identities,
the proposed modulus-nine integrality bridge, the characteristic-two quartic
identity, and the conjectured odd-index classification through degree 100.
The odd indices in this range are 0, 1, 5, 21 and 85. Successful output ends with
`ALL_CHECKS_OK`. These finite checks do not prove the unbounded conjecture.

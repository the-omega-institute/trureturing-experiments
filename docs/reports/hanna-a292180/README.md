# Hanna A292180: exact coefficient checks

The program checks Paul D. Hanna's OEIS A292180, revision 15, whose third COMMENTS line conjectures that a((2n−1)²)/4 is odd for every integer n ≥ 1. Source: <https://oeis.org/A292180/internal>; revision timestamp 2017-10-05T05:07:34-04:00; offset 1.

Run with Python 3 and the standard library:

```sh
python3 docs/reports/hanna-a292180/check.py
```

All evidence uses integers. One algorithm directly expands the source FORMULA using binomial coefficients; the other convolves powers of (1+t)/(1−t) and reflects t to −t. They agree for coefficients 0 through 1024. The program reproduces every supplied b-file value, the first 40 displayed coefficients, and all 31 square quotients printed in the source. It checks the proposed coefficient congruence on all 7,262 positive pairs (d,k) with dk ≤ 1024 and the stronger odd-index parity characterization on all 512 odd indices through 1024. The final line is ALL_CHECKS_OK. These checks establish finite agreement only; they do not prove the all-index conjecture.

Files:

- `check.py`: independent exact-arithmetic check program; SHA-256 5e7a97d5aeb807a6b9a9cc56c3f162f58b91eab7c8e0a9d66b14c6883367868d.
- `output.txt`: complete standard output from the command above.
- `b292180.txt`: unmodified OEIS source data, <https://oeis.org/A292180/b292180.txt>; SHA-256 a26d6a3a054f651f2a9c08c15b5651a66e6d9c8c279cc061186fd8d9fc70aa32.

The b-file and reproduced OEIS sequence data are attributed to Paul D. Hanna and The OEIS Foundation Inc. and retain the OEIS CC BY-SA 4.0 license: <https://oeis.org/wiki/The_OEIS_End-User_License_Agreement> and <https://creativecommons.org/licenses/by-sa/4.0/>. The check program and this README are available under CC BY-SA 4.0.

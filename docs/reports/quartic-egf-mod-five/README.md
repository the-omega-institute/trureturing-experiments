# Quartic EGF modulo five

`check.py` independently computes the exponential coefficients of the unique
zero-constant formal solution of

    A(x) = x exp(A(A(A(A(x)))))

using integer Bell recurrences through n=201. A separate ordinary-series
Fraction/Horner computation checks n=0..17. All 201 integers match the pinned
OEIS b-file, and the displayed fourth-iterate coefficients match the source.

Run with Python 3 (standard library only; network access to OEIS is required):

```sh
python3 check.py
```

The program also verifies the proposed fixed-point polynomial exactly in
F_5[X,E]/(X^5,E^5-1), compares its four iterates and exponential with independent
Hurwitz coefficient calculations through n=201, and checks the proposed
period-20 coefficient formula through n=1000. The expected final line is
`ALL_CHECKS_OK`; `output.txt` contains the successful output.

These are exact finite checks. They do not establish the all-index conjecture
in Lean. The proof still needs the integral divided-power/substitution bridge
and fixed-point identification. In particular, ordinary nilpotence alone
does not justify truncating a divided-power exponential.

Source: Paul D. Hanna, [OEIS A396804](https://oeis.org/A396804/internal),
revision #23, with composition convention and modulo-five conjecture;
[b-file](https://oeis.org/A396804/b396804.txt), SHA-256
`32a0947573161ef08b5eb7b11c581d6c676bbcba1001a39a5634c599e55ac986`.
The entry and its tabulated data are attributed to Hanna and OEIS and remain
subject to the [OEIS license agreements](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement).
The check program is new code under this repository's Apache-2.0 license.

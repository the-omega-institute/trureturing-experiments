# Odd-modulus zero classification for OEIS A397683

The conjecture says that for a positive odd integer m, the number f(m) of
consecutive units x, x+1 of equal multiplicative order modulo m is zero exactly
when m is a positive power of 3 or of 7. The source sets f(1)=1.

Source: [OEIS A397683](https://oeis.org/A397683/internal), revision 38,
by Donovan Sneider González May. The conjecture is the final Comments entry;
the prime-power classification is the first Formula entry.

Run with Python 3 and no additional dependencies:

```sh
python3 check.py
```

The stored output is `check.stdout`. All arithmetic is exact integer arithmetic.
Orders are obtained by dividing Euler's exponent and checking modular powers;
repeated multiplication independently cross-checks every unit for odd moduli
through 161. Fidelity checks reproduce the 80 published sequence terms and
the source's explicit example f(11)=4.

Full counts test the conjecture for every positive odd m through 2001 and
reproduce the source's prime-power classification for every odd prime power
through 5000. Further checks cover the exceptional local order formulas at
exponents 1 through 8 and joint primitive lifting for primes below 300 at
exponents 2 through 5.

For all positive odd m below 20000, the program separately verifies f(1)=1,
computes the zero counts at the 14 pure 3-power or 7-power exceptions, and
constructs consecutive units of maximal order at every other modulus. These
are finite exact checks; they do not establish the unbounded conjecture or
provide Lean proofs. The proposed general route uses Cohen's theorem on
consecutive primitive elements of finite fields, prime-power lifting, and CRT.

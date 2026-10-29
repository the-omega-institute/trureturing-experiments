# Lang's complex Fibonacci floor formula

Source: [OEIS A052952](https://oeis.org/A052952), revision 166. The FORMULA
conjecture is attributed to Wolfdieter Lang. The NAME is
`a(n) = Fibonacci(n+2) - (1-(-1)^n)/2`, with offset zero. The first 39 constants
in the checker are the displayed sequence; entries 39 and 40 agree with the
[linked b-file](https://oeis.org/A052952/b052952.txt).

Run with Python 3 and its standard library:

```sh
python3 docs/reports/lang-a052952/check.py
```

`output.txt` is the program's output. Every assertion passes, and its last
line is `ALL_CHECKS_OK`. The calculation uses integers and `fractions.Fraction`
throughout, with outward rounding to rational multiples of 2^-384.

The fidelity check compares the first 41 terms of the NAME with coefficients
of the stated generating function `1/((1-x-x^2)*(1-x^2))` and with the stated
inhomogeneous recurrence. This part uses integer arithmetic alone.

The settlement check independently encloses the squared complex modulus at
each integer index from 0 through 40. Put N=n+2, D=sqrt(5), phi=(1+D)/2,
x=phi^N/D, delta=1/(D*phi^N), and q=exp(-pi*N). Direct expansion of the source
definition gives

    |F(N+i*N)|^2 = x^2 + q^2*delta^2
                    - 2*(-1)^N*q*cos(2*N*log(phi))/5.

The checker uses adjacent rational bounds for sqrt(5), obtained by integer
square root; Machin's identity with alternating arctangent remainders for pi;
the positive atanh series with a geometric tail for log(phi); the exponential
series at pi/16 with a geometric tail, followed by powering and inversion;
and the cosine Taylor polynomial with its Lagrange remainder. It never uses
floating-point approximations or Binet's formula to evaluate this expression.
For the expected nonnegative integer a it verifies the rational inequalities
`a^2 <= lower_bound <= |F|^2 <= upper_bound < (a+1)^2`.

Each output row gives the integer floor and a compact rational enclosure of
`|F|^2-a^2`; the assertions use the full internal enclosures. This is a finite
independent check supporting a preregistration, not a proof for every index
and not a Lean proof. No additional package, network request, or data file is
needed to rerun it.

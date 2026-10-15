"""Actual arithmetic-symbol tail moments using existing FLINT suppliers.

Only the actual-symbol assembly, exact phase aliases and error propagation
are added here. FLINT supplies series, Bernoulli, Lerch and Hurwitz evaluation;
Johansson Eq.28 / Theorem3.3 supplies the digamma remainder. The output is
normalized scalar input, not a coupling-matrix or positivity certificate.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import time
import zlib

import flint
from flint import acb, acb_series, arb, ctx, fmpq

from arb_arithmetic import FORMULA_SOURCE, grid_interval, prime_powers


REMAINDER_SOURCE = 'https://doi.org/10.5206/mt.v3i1.14591'


def phase_representative(value, c):
    """Exact u modulo c**Z and conjugation; return representative, Im sign."""
    value = Fraction(value)
    if value <= 0 or c < 2:
        raise ValueError('Positive rational phase and c>=2 required.')
    while value < 1:
        value *= c
    while value >= c:
        value /= c
    if value == 1:
        return Fraction(1), 0
    conjugate = Fraction(c)/value
    return (value, 1) if value <= conjugate else (conjugate, -1)


class TailMoments:
    def __init__(self, c=27, M=32768, order=32):
        if c < 3 or M < 2 or order < 1:
            raise ValueError('Use integer c>=3, M>=2 and order>=1.')
        self.c, self.M, self.A, self.order = c, M, M+1, order
        self.L, self.pi = arb(c).log(), arb.pi()
        self.beta, self.a = self.pi/self.L, arb(fmpq(1, 4))
        if not self.beta > self.a:
            raise ValueError('Taylor bound requires pi/log(c)>1/4.')
        self.alpha = self.a/self.beta
        self.K = self.L*((self.L/2).cosh()-1)/self.pi
        self.delta, self.q = 1/arb(c).sqrt(), arb(fmpq(1, c*c))
        self.powers = prime_powers(c)
        weights = [arb(p).log()/arb(r).sqrt() for r, p in self.powers]
        # Reuse the existing actual-symbol uniform bound, without reoptimizing it.
        self.budget = 2*(self.L/2).cosh()+sum(weights, arb(0))
        self.primes = {}
        for (r, _), weight in zip(self.powers, weights):
            u, sign = phase_representative(r, c)
            self.primes[u] = self.primes.get(u, arb(0))+sign*weight
        self.primes = {u: w for u, w in self.primes.items() if w != 0}
        # Ordered square of the same folded P_m, including negative sum phases.
        self.square_phases = {}
        for r, vr in self.primes.items():
            for t, vt in self.primes.items():
                for phase, coefficient in [(r/t, vr*vt/2), (r*t, -vr*vt/2)]:
                    u, _ = phase_representative(phase, c)
                    self.square_phases[u] = self.square_phases.get(u, arb(0))+coefficient
        self._kernels, self._phases = {}, {}
        self._make_polynomial()

    def _make_polynomial(self):
        J, M, beta, a = self.order, self.M, self.beta, self.a
        bernoulli = [arb(fmpq.bernoulli(2*k)) for k in range(J+1)]
        # These finite rational moments are existing Lerch evaluations.
        mu = [acb(self.q).lerch_phi(-2*k, a).real for k in range(J+1)]
        old_cap = ctx.cap
        try:
            ctx.cap = 2*J
            x = acb_series([0, arb(1)/M])
            ib = acb(0, beta)
            h = (1+a*x/ib).log()-x/(2*(ib+a*x))
            ratio = x/(ib+a*x)
            for k in range(1, J):
                h -= bernoulli[k]/(2*k)*ratio**(2*k)
            # Coefficients already use t=M/m: no huge raw inverse-mode powers.
            coefficients = [-h[j].imag/2 for j in range(2*J)]
        finally:
            ctx.cap = old_cap
        coefficients[0] -= self.pi/4
        for k in range(J):
            coefficients[2*k+1] += (
                -self.K*(-self.alpha**2)**k/M**(2*k+1)
                +self.delta/2*(-1)**k*mu[k]/(beta*M)**(2*k+1))
        self.coefficients = coefficients
        self.nonzero = [(j, v) for j, v in enumerate(coefficients) if v != 0]
        self.square_coefficients = [arb(0) for _ in range(4*J-1)]
        for j, dj in self.nonzero:
            for k, dk in self.nonzero:
                self.square_coefficients[j+k] += dj*dk

        # Cauchy bound on H_J on |x|=1, before the t=M/m normalization.
        V = -(1-a/beta).log()+1/(2*(beta-a))
        for k in range(1, J):
            V += abs(bernoulli[k])/(2*k*(beta-a)**(2*k))
        # A_J and B_J divided by M**(2J) and M**(2J+1), respectively.
        self.err_a = ((1+(self.pi*(arb(J)+arb(fmpq(1, 2)))).sqrt())
                      *abs(bernoulli[J])/(4*J*(beta*M)**(2*J))
                      +V/(2*(1-arb(1)/self.A)*M**(2*J)))
        self.err_b = (self.K*self.alpha**(2*J)/M**(2*J+1)
                      +self.delta*mu[J]/(2*(beta*M)**(2*J+1)))
        self.uniform_error = self.symbol_error(self.A)

    def symbol_error(self, m):
        if m < self.A:
            raise ValueError('The common remainder applies only beyond M.')
        t = arb(self.M)/m
        return self.err_a*t**(2*self.order)+self.err_b*t**(2*self.order+1)

    def phase(self, u):
        if u not in self._phases:
            ratio = arb(fmpq(u.numerator, u.denominator))
            self._phases[u] = acb(0, 2*self.pi*ratio.log()/self.L).exp()
        return self._phases[u]

    def kernel(self, p, u=Fraction(1)):
        """M**(p-1) sum_{m>M} exp(2*pi*i*m*log_c(u))/m**p."""
        if not isinstance(p, int) or p < 2:
            raise ValueError('Integer moment power p>=2 required.')
        u, sign = phase_representative(u, self.c)
        key = (p, u)
        if key not in self._kernels:
            if u == 1:
                value = acb(p).zeta(self.A)
            else:
                z = self.phase(u)
                value = z**self.A*z.lerch_phi(p, self.A)
            value *= self.M**(p-1)
            if not value.is_finite():
                raise ArithmeticError('Nonfinite tail kernel; raise working precision.')
            self._kernels[key] = value
        value = self._kernels[key]
        return value if sign >= 0 else value.conjugate()

    def approximate_symbol(self, m):
        if m < self.A:
            raise ValueError('Use a mode beyond M.')
        t = arb(self.M)/m
        value = sum((v*t**j for j, v in self.nonzero), arb(0))
        return value-sum((w*(self.phase(u)**m).imag
                          for u, w in self.primes.items()), arb(0))

    def moments(self, p):
        """Enclose the normalized actual T_p and U_p, including all errors."""
        T = sum((v*self.kernel(p+j).real for j, v in self.nonzero), arb(0))
        T -= sum((w*self.kernel(p, u).imag for u, w in self.primes.items()), arb(0))
        U = sum((v*self.kernel(p+j).real for j, v in enumerate(self.square_coefficients)
                 if v != 0), arb(0))
        U -= 2*sum((v*w*self.kernel(p+j, u).imag
                    for j, v in self.nonzero for u, w in self.primes.items()), arb(0))
        U += sum((v*self.kernel(p, u).real for u, v in self.square_phases.items()), arb(0))
        error_T = (self.err_a*self.kernel(p+2*self.order).real
                   +self.err_b*self.kernel(p+2*self.order+1).real)
        error_U = (2*self.budget+self.uniform_error)*error_T
        return (T+arb(0, error_T.upper()), U+arb(0, error_U.upper()),
                error_T, error_U)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--c', type=int, default=27)
    parser.add_argument('--cutoff', type=int, default=32768)
    parser.add_argument('--order', type=int, default=32)
    parser.add_argument('--power', type=int, action='append')
    parser.add_argument('--precision', type=int, default=2304)
    parser.add_argument('--output-bits', type=int, default=640)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    powers = sorted(set(args.power if args.power is not None else range(2, 17)))
    if args.output_bits < 64 or args.precision < args.output_bits+64:
        parser.error('Use output-bits>=64 and at least 64 working guard bits.')
    if not powers or min(powers) < 2:
        parser.error('Request integer powers p>=2.')
    if args.out.exists():
        parser.error('Use a fresh output path.')
    ctx.prec = args.precision
    started = time.monotonic()
    arithmetic = TailMoments(args.c, args.cutoff, args.order)
    rows = []
    for p in powers:
        T, U, err_T, err_U = arithmetic.moments(p)
        rows.append({'p': p, 'T': grid_interval(T, args.output_bits),
                     'U': grid_interval(U, args.output_bits),
                     'analytic_T': grid_interval(err_T, args.output_bits+192),
                     'analytic_U': grid_interval(err_U, args.output_bits+192)})
        print('enclosed power', p, 'kernel count', len(arithmetic._kernels), flush=True)
    result = {
        'schema': 'ACTUAL_WEIL_NORMALIZED_TAIL_MOMENTS_V1',
        'c': args.c, 'M': args.cutoff, 'first_mode': args.cutoff+1,
        'definitions': {'T': 'M^(p-1)*sum_{m>M} s_m/m^p',
                        'U': 'M^(p-1)*sum_{m>M} s_m^2/m^p'},
        'formula_source': FORMULA_SOURCE, 'remainder_source': REMAINDER_SOURCE,
        'remainder_locator': 'Eq28, Theorem3.3/Eq35, derivative order1',
        'working_bits': args.precision, 'grid_bits': args.output_bits,
        'analytic_grid_bits': args.output_bits+192, 'stirling_order': args.order,
        'polynomial_degree': 2*args.order-1,
        'prime_powers': arithmetic.powers,
        'folded_sine_phases': len(arithmetic.primes),
        'folded_square_phases': len(arithmetic.square_phases),
        'normalized_polynomial': [grid_interval(v, args.output_bits)
                                  for v in arithmetic.coefficients],
        'uniform_symbol_error': grid_interval(arithmetic.uniform_error,
                                              args.output_bits+192),
        'row_count': len(rows), 'rows': rows,
        'kernel_count': len(arithmetic._kernels),
        'python': platform.python_version(), 'python_flint': flint.__version__,
        'FLINT': flint.__FLINT_VERSION__,
        'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'elapsed_seconds': time.monotonic()-started,
        'scope': 'Normalized scalar moments of the complete actual symbol, '
                 'conditional on its existing Fourier/form identification. '
                 'Not a matrix-error, retained sign, full coupling, Lean or RH certificate.'}
    raw = (json.dumps(result, sort_keys=True, separators=(',', ':'))+'\n').encode()
    args.out.write_bytes(zlib.compress(raw, level=9))
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ['rows', 'normalized_polynomial']}, indent=2))
    print('packet_sha256', hashlib.sha256(args.out.read_bytes()).hexdigest())


if __name__ == '__main__':
    main()

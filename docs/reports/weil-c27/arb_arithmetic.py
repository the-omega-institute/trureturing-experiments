"""Precision-parametric scalar inputs for the actual even Weil Fourier matrix.

FLINT supplies digamma/polygamma and directed ball arithmetic. The actual
pole-Gamma-prime formulas are reused from the project's arithmetic assembly;
this interface neither changes the form nor runs its old sign certificates.
"""
import argparse
import hashlib
import json
from math import isqrt
from pathlib import Path
import platform
import time
import zlib

import flint
from flint import acb, arb, ctx, fmpq


FORMULA_SOURCE = (
    'https://github.com/the-omega-institute/trureturing/blob/'
    '4804c4020d7d9f165a2e2683c42b01e3d47a8be4/'
    'research/weil_ground_mode/certify_prime3_refined.py')


def prime_powers(c):
    result = []
    for p in range(2, c):
        if any(p % k == 0 for k in range(2, isqrt(p)+1)):
            continue
        n = p
        while n < c:
            result.append((n, p))
            n *= p
    return sorted(result)


def grid_interval(value, bits):
    """Outward exact integer endpoints on the common 2**(-bits) grid."""
    if not value.is_finite():
        raise ArithmeticError('Nonfinite arithmetic scalar.')
    scale = 1 << bits
    lo = int((value.lower()*scale).floor().unique_fmpz())
    hi = int((value.upper()*scale).ceil().unique_fmpz())
    if hi-lo > 2:
        raise ArithmeticError('Increase working precision or series cutoff.')
    return [str(lo), str(hi)]


class Arithmetic:
    def __init__(self, c, target_bits):
        if c < 3 or target_bits < 64:
            raise ValueError('Use integer c>=3 and target_bits>=64.')
        self.c, self.d = c, arb(c).log()
        self.pi = arb.pi()
        self.H = (self.d/2).cosh()
        self.powers = prime_powers(c)
        self.primes = [(arb(p).log()/arb(n).sqrt(), arb(n).log())
                       for n, p in self.powers]
        self.q = arb(fmpq(1, c*c))
        # For b_j=2j+1/2, exp(-b_j*d)=c**(-1/2)*q**j.
        self.cutoff = 0
        threshold = arb(2)**(-target_bits-24)
        while True:
            b = arb(fmpq(4*self.cutoff+1, 2))
            decay = self.q**self.cutoff/arb(c).sqrt()
            # omega/(b*b+omega*omega)<=1/(2*b), uniformly omega>=0.
            self.symbol_tail = decay/(2*b*(1-self.q))
            # |(b-i*omega)**(-2)|<=1/b**2.
            self.diagonal_tail = 2*decay/(self.d*b*b*(1-self.q))
            if self.symbol_tail < threshold and self.diagonal_tail < threshold:
                break
            self.cutoff += 1
        self.gamma_terms = []
        decay = 1/arb(c).sqrt()
        for j in range(self.cutoff):
            self.gamma_terms.append((arb(fmpq(4*j+1, 2)), decay))
            decay *= self.q

    def scalars(self, n):
        if n < 0:
            s, a = self.scalars(-n)
            return -s, a
        omega = 2*self.pi*n/self.d
        z = acb(arb(fmpq(1, 4)), omega/2)
        psi, tri = z.digamma(), z.polygamma(1)
        correction_s, correction_a = arb(0), arb(0)
        for b, decay in self.gamma_terms:
            correction_s += omega*decay/(b*b+omega*omega)
            correction_a += decay*(1/acb(b, -omega)**2).real
        symbol = (-2*omega*(self.H-1)/(omega*omega+arb(fmpq(1, 4)))
                  -psi.imag/2+correction_s)
        symbol += arb(0, self.symbol_tail.upper())
        diagonal = (psi.real-self.pi.log()+tri.real/(2*self.d)
                    -2*correction_a/self.d
                    +4*(self.H-1)/self.d
                    *(1/acb(arb(fmpq(1, 2)), omega)**2).real)
        diagonal += arb(0, self.diagonal_tail.upper())
        for weight, shift in self.primes:
            symbol -= weight*(omega*shift).sin()
            diagonal -= 2*weight*(1-shift/self.d)*(omega*shift).cos()
        # The exact zero symbol is zero; do not attach a spurious tail there.
        return (arb(0) if n == 0 else symbol), diagonal

    def even_entry(self, m, n, scalars):
        """Orthonormal even basis: 1/sqrt(d), sqrt(2/d)*cos(2*pi*n*x/d)."""
        if m < 0 or n < 0:
            raise ValueError('Nonnegative even mode indices required.')
        if m < n:
            return self.even_entry(n, m, scalars)
        sm, am = scalars[m]
        sn, _ = scalars[n]
        if m == n:
            return am if n == 0 else am-sn/(self.pi*n)
        phase = -1 if (m+n) % 2 else 1
        if n == 0:
            return -phase*arb(2).sqrt()*sm/(self.pi*m)
        return phase*2*(n*sn-m*sm)/(self.pi*(m*m-n*n))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--c', type=int, default=27)
    parser.add_argument('--last-low', type=int, default=1280)
    parser.add_argument('--extra-mode', type=int, action='append', default=[])
    parser.add_argument('--precision', type=int, default=768)
    parser.add_argument('--output-bits', type=int, default=640)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.precision < args.output_bits+64 or args.last_low < 0:
        parser.error('Use at least 64 working guard bits and nonnegative modes.')
    if any(n < 0 for n in args.extra_mode):
        parser.error('Extra modes must be nonnegative.')
    if args.out.exists():
        parser.error('Use a fresh output path.')
    ctx.prec = args.precision
    started = time.monotonic()
    arithmetic = Arithmetic(args.c, args.output_bits)
    rows = []
    for n in sorted(set(range(args.last_low+1)) | set(args.extra_mode)):
        s, a = arithmetic.scalars(n)
        rows.append({'n': n, 's': grid_interval(s, args.output_bits),
                     'a': grid_interval(a, args.output_bits)})
    result = {
        'schema': 'ACTUAL_WEIL_ARITHMETIC_SCALARS_V1',
        'c': args.c, 'full_length': 'log(c)', 'last_low': args.last_low,
        'basis': '1/sqrt(d),sqrt(2/d)*cos(2*pi*n*x/d); no (-1)^n phase',
        'formula_source': FORMULA_SOURCE,
        'prime_powers': arithmetic.powers,
        'working_bits': args.precision, 'grid_bits': args.output_bits,
        'gamma_correction_terms': arithmetic.cutoff,
        'symbol_tail_upper': grid_interval(arithmetic.symbol_tail, args.output_bits+32),
        'diagonal_tail_upper': grid_interval(arithmetic.diagonal_tail, args.output_bits+32),
        'tail_grid_bits': args.output_bits+32,
        'row_count': len(rows), 'rows': rows,
        'python': platform.python_version(), 'python_flint': flint.__version__,
        'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'elapsed_seconds': time.monotonic()-started,
        'scope': 'Actual scalar assembly only, conditional on the existing '
                 'Fourier/form identification. No low-block sign, infinite '
                 'coupling, positivity, Lean or RH certificate.'}
    raw = (json.dumps(result, sort_keys=True, separators=(',', ':'))+'\n').encode()
    args.out.write_bytes(zlib.compress(raw, level=9))
    print(json.dumps({k: v for k, v in result.items() if k != 'rows'}, indent=2))
    print('packet_sha256', hashlib.sha256(args.out.read_bytes()).hexdigest())


if __name__ == '__main__':
    main()

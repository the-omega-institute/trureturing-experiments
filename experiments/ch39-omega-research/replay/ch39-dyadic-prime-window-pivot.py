"""Enclose three-node Schur readings of one or two dyadic prime windows."""

from flint import acb, arb, arb_series, ctx


def arithmetic_times():
    half = arb(1) / 2
    s = arb_series([half, 1], prec=3)
    xi = s * (s - 1) / 2 * (-s * arb.pi().log() / 2).exp() * (s / 2).gamma() * s.zeta()
    mass = xi[2] / xi[0]
    readings = []
    for j in (1, 2):
        s = half + j
        z = arb_series([s, 1], prec=2).zeta()
        log_deriv = 1/s + 1/(s-1) + ((s/2).digamma()-arb.pi().log())/2 + z[1]/z[0]
        readings.append(log_deriv / (2*j*mass))
    return readings


def zero_data(start, count, log_two, root_two):
    rows = []
    for zero in acb.zeta_zeros(start, count):
        g = zero.imag
        x = 1 / g**2
        bare = 1 / (g**2 + arb(1)/4)
        cosine = (g * log_two).cos()
        residue = bare * (3*root_two - 4*cosine)
        first = x / (1+x)
        second = 4*x*x / ((1+x)*(1+4*x))
        rows.append((cosine, bare, residue, first, second))
    return rows


def weighted_sums(rows, k, root_two):
    a = b = c = arb(0)
    constant = 1 + 2*k*k
    slope = 2*root_two*k
    for cosine, _, residue, first, second in rows:
        weight = residue * (constant + slope*cosine).sqrt()
        a += weight
        b += weight * first
        c += weight * second
    return a, b, c


def schur_readout(q, z):
    # q=B/A, z=C/B for the first and second differences of the time vectors.
    residual = q*(1-z) * (z - 4*q*(1-z)/(1-q))
    assert residual > 0
    return ((1-4*q+3*q*z)/residual).sqrt()


def verify():
    ctx.prec, ctx.cap = 120, 3
    root_two = arb(2).sqrt()
    log_two = arb(2).log()
    rows = zero_data(1, 1024, log_two, root_two)
    six = zero_data(11485, 6, log_two, root_two)
    total = 1 + arb.const_euler()/2 - (4*arb.pi()).log()/2
    tail = total - sum(row[1] for row in rows)
    assert tail > 0
    first_max, second_max = rows[-1][-2:]
    t1, t2 = arithmetic_times()
    q_target = 1 - t1

    # For every positive-zero probability measure, z is a weighted average
    # of 4/(gamma^2+4). At fixed q the threshold polynomial decreases here.
    upper_z = arb(1)/50
    assert 4/(acb.zeta_zero(1).imag**2+4) < upper_z
    target = arb('266.6444299')
    def threshold_gap(z):
        residual = q_target*(1-z)*(z-4*q_target*(1-z)/(1-q_target))
        return 1-4*q_target+3*q_target*z-target**2*residual

    derivative_upper = 3-target**2*(1-2*upper_z)
    assert derivative_upper < 0
    low, high = arb(0), upper_z
    assert threshold_gap(low) > 0 > threshold_gap(high)
    for _ in range(72):
        middle = (low+high)/2
        sign = threshold_gap(middle)
        if sign > 0:
            low = middle
        else:
            assert sign < 0
            high = middle
    crossing = 1-4*q_target+3*q_target*low.union(high)
    assert arb('2.16095648e-7') < t2-crossing < arb('2.16095649e-7')
    print('T2 boundary at the arithmetic T1:', crossing.str(20))

    def moments(k):
        a, b, c = weighted_sums(rows, k, root_two)
        upper = (3*root_two+4)*(1+abs(k)*root_two)*tail
        lower = (3*root_two-4)*abs(1-abs(k)*root_two)*tail
        return (a + lower.union(upper), b + arb(0).union(upper*first_max),
                c + arb(0).union(upper*second_max))

    for k, pivot_lo, pivot_hi, ppm_lo, ppm_hi in (
        (0, 239, 251, '2.22', '2.34'), (1, 273, 300, '2.07', '2.31')
    ):
        a, b, c = moments(arb(k))
        # Keep A explicit here to avoid using its tail enclosure multiple times.
        residual = c*(1-c/b) - 4*(b-c)**2/(a-b)
        assert residual > 0
        readout = ((a-4*b+3*c)/residual).sqrt()
        ppm = 10**6 * weighted_sums(six, arb(k), root_two)[0]/a
        assert pivot_lo < readout < pivot_hi
        assert arb(ppm_lo) < ppm < arb(ppm_hi)
        print(f'k={k}: pivot', readout.lower().str(12), readout.upper().str(12),
              'band in ppm', ppm.lower().str(9), ppm.upper().str(9))

    left, right = arb(34)/100, arb(1)/2
    for k, sign in ((left, -1), (right, 1)):
        a, b, _ = moments(k)
        assert sign*(q_target - b/a) > 0

    # At an exact T1 match q is known. Bound C/B on every parameter cell.
    all_pivots = all_gaps = None
    for j in range(320):
        lo = left + j/arb(2000)
        k = lo.union(lo + 1/arb(2000))
        _, b, c = moments(k)
        z = c/b
        readout = schur_readout(q_target, z)
        gap = 1-4*q_target+3*q_target*z-t2
        band_ppm = 10**6 * weighted_sums(six, k, root_two)[0] * q_target/b
        assert 0 < band_ppm < 3
        assert arb('266.8') < readout < 285
        assert arb('-5.52e-6') < gap < arb('-2.74e-7')
        all_pivots = readout if all_pivots is None else all_pivots.union(readout)
        all_gaps = gap if all_gaps is None else all_gaps.union(gap)
    print('Enclosure at exact T1 match:')
    for name, value in [('pivot', all_pivots), ('T2 minus RH', all_gaps)]:
        print(name, value.lower().str(12), value.upper().str(12))


if __name__ == '__main__':
    verify()

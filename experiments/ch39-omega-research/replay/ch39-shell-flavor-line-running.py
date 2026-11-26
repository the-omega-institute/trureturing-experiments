"""One-loop motion of the flavor lines required by shell exchange protection."""
from pathlib import Path
from runpy import run_path
from flint import acb_mat, arb, ctx

mixing = run_path(str(Path(__file__).with_name('ch39-flavour-symmetry-gap.py')))


def squares(values):
    return acb_mat([[v*v if i == j else 0 for j in range(3)] for i, v in enumerate(values)])


def angular_speed(opposite, singular, j, right):
    terms = []
    for i in range(3):
        if i == j:
            continue
        weight = singular[i]*singular[j] if right else (singular[i]**2 + singular[j]**2)/2
        terms.append((abs(opposite[i, j])*weight/(singular[j]**2 - singular[i]**2))**2)
    return 3*sum(terms).sqrt()/(16*arb.pi()**2)


def main():
    ctx.prec = 192
    v = acb_mat(mixing['ckm'](*[arb(a).sin() for a in ('.22704', '.04275', '.00377')],
                              arb('1.139')))
    up = list(map(arb, ('.00000615', '.00311', '.8616')))
    down = list(map(arb, ('.0000135', '.000268', '.01401')))
    star = v.conjugate().transpose()
    up_in_down = star*squares(up)*v
    down_in_up = v*squares(down)*star
    comm = up_in_down*squares(down) - squares(down)*up_in_down
    cp = abs(comm.det())
    assert arb('9.43038e-25') < cp < arb('9.43039e-25')
    dr = angular_speed(up_in_down, down, 0, True)
    ur = angular_speed(down_in_up, up, 0, True)
    bottom = angular_speed(up_in_down, down, 2, False)
    assert arb('2.88958e-7') < dr < arb('2.88960e-7')
    assert arb('1.53481e-12') < ur < arb('1.53483e-12')
    assert arb('0.00030245858') < bottom - ur < bottom + ur < arb('0.00030245859')
    cross_gram = acb_mat([[down[i]*down[j]*up_in_down[i,j] for j in range(3)]
                          for i in range(3)])
    trace = cross_gram.trace().real
    moment = cross_gram[0,0].real / trace
    derivative = 6*sum(abs(cross_gram[j,0])**2/(down[j]**2 - down[0]**2)
                       for j in (1,2))/(16*arb.pi()**2*trace)
    assert arb('7.17391e-11') < moment < arb('7.17392e-11')
    assert arb('2.00603e-12') < derivative < arb('2.00604e-12')
    print('Light right-line invariant expectation', moment.str(16),
          'increase per unit log scale', derivative.str(16))
    print('Commutator determinant modulus', cp.str(16))
    print('Per unit log scale: d singlet', dr.str(16), 'u singlet', ur.str(16),
          'left bottom', bottom.str(16))


if __name__ == '__main__':
    main()

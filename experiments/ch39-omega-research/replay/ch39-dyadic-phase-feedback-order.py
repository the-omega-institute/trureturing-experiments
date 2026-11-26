"""Sharp rational order for a dyadic clock on the six-zero shell band."""
from flint import acb, arb, arb_mat, ctx, fmpz_poly


def positive_definite(matrix):
    for size in range(1, matrix.nrows()+1):
        leading = arb_mat([[matrix[i,j] for j in range(size)] for i in range(size)])
        assert leading.det() > 0


def verify():
    ctx.prec = 160
    heights = [z.imag for z in acb.zeta_zeros(11485, 6)]
    log_two = arb(2).log()
    delta = [[x-y for y in heights[1::2]] for x in heights[::2]]
    loewner = arb_mat([[2*(d*log_two/2).sin()/d for d in row] for row in delta])
    error = arb_mat([[2/abs(d) for d in row] for row in delta])
    identity = arb_mat([[int(i==j) for j in range(3)] for i in range(3)])
    s_min, e_max = arb('0.00086'), arb('6.49')
    positive_definite(loewner.transpose()*loewner - identity*s_min**2)
    positive_definite(identity*e_max**2 - error.transpose()*error)
    tolerance = arb(1)/8192
    assert s_min/e_max > tolerance

    center = arb(11137)
    slope = arb(45426)/65536
    r = -arb(46905692)/2**24
    fixed = acb(1,r)/acb(1,-r)
    width = max(center-heights[0], heights[-1]-center)
    rotation_error = abs(acb(0,center*log_two).exp()-fixed)
    bound = rotation_error + width*abs(slope-log_two) + (width*slope)**7/100800
    assert bound < arb('0.000098') < tolerance
    for g in heights:
        z = acb(0, (g-center)*slope)
        p = 120 + z*(60 + z*(12+z))
        assert abs(fixed*p/p.conjugate() - acb(0, g*log_two).exp()) < bound

    g, h = heights[-2:]
    arc = 2*((arb(1)/(2*g)).atan() - (arb(1)/(2*h)).atan())
    stretch = (2*(((h-g)*log_two/2).sin())-2*tolerance)/arc
    critical_gap = 6/(stretch+3)
    assert critical_gap < arb('7.02e-8')
    gaps = []
    for z, multiplicity in fmpz_poly([120,60,12,1]).complex_roots():
        assert multiplicity == 1 and z.real < 0
        g = center+z/acb(0,slope)
        a = (g-acb(0,arb(1)/2))/(g+acb(0,arb(1)/2))
        gaps.append(1-abs(a))
    assert all(arb('4.27e-8') < gap < arb('5.41e-8') for gap in gaps)
    print('Degree-two uniform error lower bound:', (s_min/e_max).lower().str(12))
    print('Degree-three full-band uniform error upper bound:', bound.upper().str(12))
    print('Largest feedback radius has gap below', critical_gap.upper().str(12))
    print('Three realized radial gaps:', [x.str(12) for x in gaps])


if __name__ == '__main__':
    verify()

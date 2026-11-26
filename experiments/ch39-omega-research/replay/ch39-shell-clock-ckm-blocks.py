"""Certify the alternating phase-block count on the two shell arc normals."""
from flint import acb, acb_mat, arb, ctx

ctx.prec = 192
radial = (1, 5, 29)
centers = tuple(map(arb, (".22501", ".04183", ".003732", "1.147")))
sine_radius = arb(".000055")
phase_radius = arb(".006")


def cross(x, y):
    return tuple(x[(i+1) % 3]*y[(i+2) % 3] - x[(i+2) % 3]*y[(i+1) % 3]
                 for i in range(3))


def dot(x, y):
    return sum(a*b for a, b in zip(x, y))


def cis(x):
    s, c = x.sin_cos()
    return acb(c, s)


axes = tuple(cross(radial, endpoint) for endpoint in ((13, 23, 13), (23, 13, 13)))
fixed_lines = tuple(cross(radial, axis) for axis in axes)
axis_cos_sq = arb(dot(*axes)**2) / (dot(axes[0], axes[0])*dot(axes[1], axes[1]))


def target_box(sine_halfwidth=sine_radius, phase_halfwidth=phase_radius):
    sines = [v + arb(0, sine_halfwidth.upper()) for v in centers[:3]]
    return standard_matrix(*sines, centers[3] + arb(0, phase_halfwidth.upper()))


def standard_matrix(s, t, u, phase):
    a, b, c = ((1-v*v).sqrt() for v in (s, t, u))
    d = cis(phase)
    return ((a*c, s*c, u*d.conjugate()),
            (-s*b-a*t*u*d, a*b-s*t*u*d, t*c),
            (s*t-a*b*u*d, -a*t-s*b*u*d, b*c))


def rephasing_equations(v):
    # The first two output modulus constraints have only three interference terms.
    offsets = []
    coeffs = []
    for i, row in enumerate(v[:2]):
        scaled = [acb(entry)*weight for entry, weight in zip(row, radial)]
        offsets.append(sum(abs(z)**2 for z in scaled) - radial[i]**2)
        coeffs.append([2*scaled[a]*scaled[b].conjugate() for a, b in ((0,1), (0,2), (1,2))])
    return offsets, coeffs


def exclude_three():
    v = target_box()
    offsets, coeffs = rephasing_equations(v)
    pending = [(0, 0, 0)]
    boxes = 0
    depth = 0
    while pending:
        level, x, y = pending.pop()
        boxes += 1
        depth = max(depth, level)
        unit = arb.pi() / (1 << level)
        phi, psi = ((2*i+1)*unit + arb(0, unit.upper()) for i in (x,y))
        h, k = cis(phi), cis(psi)
        diffs = (cis(phi-psi), h, k)
        if any(0 not in offset + sum((a*b).real for a,b in zip(row,diffs))
               for offset, row in zip(offsets, coeffs)):
            continue

        right = (h, k, acb(1))
        weighted = [tuple(z*phase for z, phase in zip(row, right)) for row in v]
        left = [(dot(row, radial)/r).conjugate() for row,r in zip(weighted, radial)]
        possible = False
        for m in fixed_lines:
            value = sum(mi*li*dot(row, m) for mi,li,row in zip(m,left,weighted))/dot(m,m)
            if 0 in abs(value - axis_cos_sq) - (1 - axis_cos_sq):
                possible = True
                break
        if possible:
            assert level < 19, (level,x,y)
            pending.extend((level+1, 2*x+a, 2*y+b) for a in (0,1) for b in (0,1))
    return boxes, depth


angles = [(1-arb(d)/1734).acos() for d in (244,724,804)]
counts = (4227019212963, 4682583908189, 3807923366500, 3273322617774)


def clock_word(height):
    clocks = [height*(angle-angles[0]) for angle in angles[1:]]
    identity = acb_mat([[int(i == j) for j in range(3)] for i in range(3)])
    v = identity
    for axis, angle, count in zip(axes*2, clocks*2, counts):
        factor = (cis(count*angle)-1) / dot(axis,axis)
        v = v * (identity + acb_mat([[x*y*factor for y in axis] for x in axis]))
    return v


def check_match(index, v):
    u = abs(v[0,2])
    c = (1-u*u).sqrt()
    s, t = abs(v[0,1])/c, abs(v[1,2])/c
    a, b = (1-s*s).sqrt(), (1-t*t).sqrt()
    j = (v[0,0]*v[1,0].conjugate()*v[1,1]*v[0,1].conjugate()).imag
    sin_delta = j / (a*b*c*c*s*t*u)
    cos_delta = (abs(v[1,0])**2 - (s*b)**2 - (a*t*u)**2) / (2*a*b*s*t*u)
    assert sin_delta > 0 and cos_delta > 0
    delta = (sin_delta/cos_delta).atan()
    errors = [abs(x-y) for x,y in zip((s,t,u,delta),centers)]
    assert all(e < arb(".00000247") for e in errors[:3])
    assert errors[3] < arb(".000227")
    print(index, 'errors', *(e.upper().str(6) for e in errors), 'J', j.str(13))


if __name__ == "__main__":
    print('rephasing boxes and maximum depth', *exclude_three())
    matches = []
    for index, zero in enumerate(acb.zeta_zeros(10815, 1278), 10815):
        v = clock_word(zero.imag)
        u = abs(v[0,2])
        if abs(u - centers[2]) > sine_radius:
            continue
        c = (1-u*u).sqrt()
        if any(abs(vij/c-center) > sine_radius
               for vij, center in zip((abs(v[0,1]), abs(v[1,2])), centers)):
            continue
        check_match(index, v)
        matches.append(index)
    assert matches == [11439, 11737]
    print('matching indices in 10815..12092:', *matches)

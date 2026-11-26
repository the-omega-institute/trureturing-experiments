"""三低模、单位面强度下四维单边 Haar 刚度的无截频求积。"""
from math import lcm

import numpy as np
from flint import arb, ctx


def tables(period):
    face, slope, weyl = (np.empty(period) for _ in range(3))
    for i in range(period):
        if i in (0, period // 2):
            values = (1, 0, 0) if i == 0 else (0, 0, 4)
        else:
            sine, cosine = (arb.pi() * i / period).sin_cos()
            values = (cosine**36, sine/cosine, 4*sine*sine)
        for arr, v in zip((face, slope, weyl), values):
            arr[i] = float(v)
            if isinstance(v, arb):
                assert abs(arb(arr[i])-v) < abs(v)*arb("1e-15")
    assert np.max(np.abs(slope)) < 1100
    return face, slope, weyl


def grid_moments(ny, nw, nc, by_charge, weyl, absolute_gradients=False):
    assert ny % 6 == nw % 2 == nc % 3 == 0
    period = lcm(ny, nw, nc)
    # All six permutations have the same squared Cartan gradient. Coincident
    # colours vanish in the Weyl density. Weak endpoints vanish as well.
    triples = [(a,b,(-a-b) % nc) for a in range(nc) for b in range(a+1,nc)
               if b < (-a-b) % nc]
    colors = (np.array(triples).T * (period//nc))[:, None, :]
    weak = (np.arange(1,nw//2) * (period//nw))[:, None]
    assert weak.size * len(triples) < 1000000
    haar = weyl[2*weak] * np.prod([weyl[(colors[i]-colors[j]) % period]
                                    for i,j in ((0,1),(1,2),(0,2))], axis=0)
    sums = np.zeros(4)
    for i in range(ny//6):
        y = i * (period//ny)
        density = haar.copy()
        grad_y, grad_w = np.zeros_like(haar), np.zeros_like(haar)
        grad_c = []
        # Three colour eigenangles. Each carries the two Q weights and the
        # conjugate u,d weights. The slope tables give the angular gradients.
        for c in colors:
            grad = np.zeros_like(haar)
            for q, w, sign in ((1,1,1),(1,-1,1),(-4,0,-1),(2,0,-1)):
                face, slope = by_charge[abs(q)]
                angle = (q*y + w*weak + sign*c) % period
                density *= face[angle]
                grad_y += (abs(q) if absolute_gradients else q)*slope[angle]
                grad_w += (abs(w) if absolute_gradients else w)*slope[angle]
                grad += (1 if absolute_gradients else sign)*slope[angle]
            grad_c.append(grad)
        for q, w in ((-3,1),(-3,-1),(6,0)):
            face, slope = by_charge[abs(q)]
            angle = (q*y + w*weak) % period
            density *= face[angle]
            grad_y += (abs(q) if absolute_gradients else q)*slope[angle]
            grad_w += (abs(w) if absolute_gradients else w)*slope[angle]
        color_sign = 1 if absolute_gradients else -1
        color_square = sum((grad_c[a]+color_sign*grad_c[b])**2 for a,b in ((0,1),(0,2),(1,2)))
        sums += [np.sum(density), np.sum(density*grad_y**2),
                 np.sum(density*grad_w**2), np.sum(density*color_square)]
    # Recover weighted RMS values. Outward allowances cover binary64 table
    # and product errors, nonnegative sums, and absolute gradient errors.
    assert np.finfo(float).eps == 2**-52 and sums[0] > 1e-9
    norm = [(arb(float(value))/arb(float(sums[0]))).sqrt() for value in sums[1:]]
    enclosed = []
    for x in norm:
        err = arb("2e-8")*(1+x)
        enclosed.append((x-err).union(x+err))
    weight = arb(float(sums[0]))*72/(ny*nw*nc**2)
    weight *= arb("1 +/- 1e-9")
    return weight, enclosed


def one_loop(ratios):
    normalized_y = ratios[0]/60
    alpha_em_inverse, weak_angle = arb("127.930"), arb(".23129")
    a1 = 3*(1-weak_angle)*alpha_em_inverse/5
    a2 = weak_angle*alpha_em_inverse
    b1, b2, b3 = arb(41)/10, -arb(19)/6, -7
    clock = (a1 - normalized_y*a2)/(b1-normalized_y*b2)
    common = a2-b2*clock
    scale = arb("91.1876")*(2*arb.pi()*clock).exp()
    return scale, 1/(ratios[1]*common+b3*clock)


def verify():
    ctx.prec = 160
    ny, nw, nc = 660, 150, 150
    # A factor cos(phi/2)^36 has Laurent degree 18. Each differentiated
    # factor retains that degree, including cross terms of the gradients.
    assert ny > 18*36 and nw > 18*8+2 and nc > 18*8+4
    face, slope, weyl = tables(lcm(ny,nw,nc))
    pairs = dict.fromkeys((1,2,3,4,6), (face, slope))
    _, (y,w,c) = grid_moments(ny, nw, nc, pairs, weyl)
    weak = 45*w*w/12
    ratios = [45*y*y/weak, (45*c*c/48)/weak]
    print("三模单位强度的共享边刚度 Y、W、C", 45*y*y, weak, 45*c*c/48)
    print("三模单位强度的共享边超荷/弱、色/弱", *ratios)
    assert arb("60.3222") < ratios[0] < arb("60.3224")
    assert arb(".993111") < ratios[1] < arb(".993113")
    # Conditional matching prescription: use this single-link Hessian as
    # the entire UV boundary, then run only the minimal SM at one loop.
    scale, strong = one_loop(ratios)
    assert arb(".07194") < strong < arb(".07195") and arb("8.32e12") < scale < arb("8.34e12")
    print("一圈匹配能标 GeV、强耦合", scale, strong)
    return ratios


if __name__ == "__main__":
    verify()

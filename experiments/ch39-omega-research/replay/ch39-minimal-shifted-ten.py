"""独立物质数连接的最低维数，以及借用原生输运的一份零裸质量十维。"""
from itertools import combinations, product
from fractions import Fraction
from math import prod
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_poly, ctx

crossing = run_path(str(Path(__file__).with_name("ch39-unit-seventy-five-crossing.py")))
soft, smooth = crossing["soft"], crossing["smooth"]


def verify():
    ctx.prec = 320
    def dimension(label):
        return prod(Fraction(j-i+sum(label[i:j]),j-i) for i in range(4) for j in range(i+1,5))
    assert all(dimension(tuple(2*int(i==j) for i in range(4))) >= 10 for j in range(4))
    assert sorted(dimension(l) for l in product(range(2),repeat=4) if dimension(l)<10) == [1,5,5]
    masses = crossing["shell_spectra"]((0,5))
    geo = soft["geo"]
    r = geo["theta"](804)/geo["theta"](12)
    mz,em,weak = arb("91.1876"),arb("127.930"),arb(".23129")
    initial = (3*(1-weak)*em/5, weak*em, 1/arb(".1180"))
    pi, beta = arb.pi(), smooth["beta"]
    flat = sum(smooth["sheet"](v,r*r,False) for v in masses[0])
    window = (arb("1e17")/mz).log()
    cmax = initial[0]-3*initial[1]+2*initial[2]+window/(5*pi)
    five = arb("2.4")*(flat+window)/(3*pi)
    assert cmax+five < -arb(".1972")
    print("Universal single five projected gap upper bound", (cmax+five).upper().str(11))

    h = acb.zeta_zero(11431).imag
    star, ordinary, _ = soft["vector_pencil"](h)
    loss = arb(23)/25*arb(2)**-20/(h*h+soft["point"]["study"]["offset"])
    pole, orbit = soft["pole"],soft["star"]["orbit"]
    def characteristic(sign):
        m = star+ordinary
        for j in orbit:
            m[pole,pole] -= loss
            m[j,j] -= loss
            m[pole,j] = m[j,pole] = -(1-loss)*sign
        return arb_poly([z.real for z in m.charpoly().coeffs()])
    plus, minus = characteristic(1),characteristic(-1)
    assert plus[0].contains(0)
    plus = arb_poly([0,*plus.coeffs()[1:]])
    def one_minus_cos(x):
        return -2*plus(x)/(minus(x)-plus(x))
    # Parameterize by the lowest transported weak-doublet/color-triplet
    # eigenvalue x (the bare mass is exactly zero).
    mass_box = arb("5.8317e-19").union(arb("5.8319e-19"))
    t = one_minus_cos(mass_box)
    assert 0 < t < arb("1e-15")
    motion = arb(0,(arb(6)*t).sqrt().upper())
    mu_box = arb("6.3573e15").union(arb("6.3575e15"))
    common,light_mass = matching_box(mass_box,motion,mu_box)
    angle = 2*(t/2).sqrt().asin()
    assert arb("2.2022e-9") < angle < arb("2.2024e-9")
    assert arb("53.9033") < common < arb("53.9035")
    assert arb("5.687e5") < light_mass < arb("5.689e5")
    print("Ten matching angle, scale, inverse, lightest physical mass GeV",
          angle.str(10),mu_box.str(10),common.str(10),light_mass.str(10))



def matching_box(mass_box,motion,mu_box):
    masses = crossing["shell_spectra"]((0,5))
    geo = soft["geo"]
    r = geo["theta"](804)/geo["theta"](12)
    mz,em,weak = arb("91.1876"),arb("127.930"),arb(".23129")
    initial = (3*(1-weak)*em/5, weak*em, 1/arb(".1180"))
    pi, beta = arb.pi(), smooth["beta"]
    # T_Y,T_W,T_C in wedge^2(5), grouped by signed hypercharge times six.
    weight = {}
    for i,j in combinations(range(5),2):
        labels = [int(k in (i,j)) for k in range(5)]
        q = -2*sum(labels[:3])+3*sum(labels[3:])
        row = weight.setdefault(q,[arb(0) for _ in range(3)])
        row[0] += arb(q*q)/60
        row[1] += arb((labels[3]-labels[4])**2)/4
        row[2] += arb((labels[0]-labels[1])**2)/4
    assert set(weight) == {-4,1,6}
    assert all(sum(row[i] for row in weight.values()).contains(arb("1.5")) for i in range(3))

    towers = {q:[v+motion for v in masses[abs(q-1)]] for q in weight}
    assert all(v > 1 for row in towers.values() for v in row)
    vect = [sum(smooth["sheet"](v,r*r,True) for v in masses[q]) for q in (0,5)]
    bg = [a+7*(n*vect[0]+b*vect[1])/(4*pi) for a,n,b in zip(initial,(0,2,3),(5,3,2))]
    def intercept(x):
        d = {q:sum(smooth["sheet"](v,r*r,False) for v in row) for q,row in towers.items()}
        # Monotonicity gives tight range even though the input x is tiny.
        d[1] += smooth["sheet"](x.lower(),r*r,False).union(smooth["sheet"](x.upper(),r*r,False))
        return [a-2*sum(row[i]*d[q] for q,row in weight.items())/(3*pi) for i,a in enumerate(bg)]
    def separate(a):
        return (beta[1]-beta[2])*(a[0]-a[2])-(beta[0]-beta[2])*(a[1]-a[2])
    logs = (mu_box/mz).log()
    ref2 = r*r*(-2*logs).exp()
    error = arb(0,(ref2/8*(19*(35/(4*pi)+1/pi)+1/(pi*mass_box.lower()))).upper())
    def read(x,L):
        return [a-b*L/(2*pi)+error for a,b in zip(intercept(x),beta)]
    faces = [separate(read(x,logs.mid())) for x in (mass_box.lower(),mass_box.upper())]
    assert faces[0]*faces[1] < 0, faces
    for L,sign in ((logs.lower(),1),(logs.upper(),-1)):
        vals = read(mass_box,L)
        assert sign*(vals[0]-vals[2]) > 0
    common = read(mass_box,logs)[2]
    def vector_lower(row):
        # For f>=1 the vector slope is >= its f=1 lower Riemann integral.
        total = arb(0)
        for j in range(32):
            lo,hi = arb(j)/32,arb(j+1)/32
            w = (8*(hi-lo)-hi**3+lo**3)/7
            z = r*r*(1-hi*hi)/4
            total += w*sum(z/(v+z) for v in row)
        return total
    neutral,charged = (vector_lower(masses[q]) for q in (0,5))
    slopes = [(-b/2+7*(n*neutral+k*charged)/4-19)/pi
              for b,n,k in zip(beta,(0,2,3),(5,3,2))]
    assert all(z>arb("17.68") for z in slopes)
    print("Post-matching uniform positive log-momentum slopes",*[z.lower().str(8) for z in slopes])

    return common,mass_box.sqrt()*mu_box/r


if __name__ == "__main__":
    verify()

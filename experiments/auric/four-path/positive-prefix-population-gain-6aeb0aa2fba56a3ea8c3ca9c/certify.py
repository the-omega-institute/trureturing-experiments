"""Exact continuous certificates for the FourPath positive-prefix gain bound.

Only Python integer arithmetic decides assertions. Every real interval endpoint
is an integer divided by S=2**96; operations round outwards. The certificate
manifest retains complete closed leaves, expression images and thresholds.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as Q
from hashlib import sha256
from math import isqrt
from pathlib import Path
import json

S = 1 << 96


def ceildiv(a, b):
    return -((-a) // b)


@dataclass(frozen=True)
class I:
    lo: int
    hi: int

    def __post_init__(self):
        assert self.lo <= self.hi

    @classmethod
    def rational(cls, value):
        f = Q(value)
        return cls(f.numerator*S//f.denominator,
                   ceildiv(f.numerator*S, f.denominator))

    @classmethod
    def box(cls, a, b):
        assert Q(a) <= Q(b)
        return cls(cls.rational(a).lo, cls.rational(b).hi)

    @staticmethod
    def coerce(other):
        return other if isinstance(other, I) else I.rational(other)

    def __add__(self, other):
        b = self.coerce(other)
        return I(self.lo+b.lo, self.hi+b.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-self.coerce(other))

    def __rsub__(self, other):
        return self.coerce(other) + (-self)

    def __mul__(self, other):
        b = self.coerce(other)
        products = [a*c for a in (self.lo, self.hi) for c in (b.lo, b.hi)]
        return I(min(products)//S, ceildiv(max(products), S))

    __rmul__ = __mul__

    def recip(self):
        assert self.lo > 0 or self.hi < 0, (self.lo, self.hi)
        return I(S*S//self.hi, ceildiv(S*S, self.lo))

    def __truediv__(self, other):
        return self*self.coerce(other).recip()

    def __rtruediv__(self, other):
        return self.coerce(other)*self.recip()

    def sq(self):
        a = 0 if self.lo <= 0 <= self.hi else min(self.lo*self.lo, self.hi*self.hi)
        b = max(self.lo*self.lo, self.hi*self.hi)
        return I(a//S, ceildiv(b, S))

    def sqrt(self):
        assert self.lo >= 0, (self.lo, self.hi)
        a = isqrt(self.lo*S)
        b = isqrt(self.hi*S)
        return I(a, b if b*b == self.hi*S else b+1)


h = I.rational('1/2').sqrt()
c = ((1+h)/2).sqrt()
s2 = (1-h)/2
s = s2.sqrt()
k = c-h
sqrt3 = I.rational(3).sqrt()


def data(xi, q):
    x = 1-xi/4+q/2
    y = (1-x.sq()).sqrt()
    u = (1-q).sqrt()
    v = q.sqrt()
    R = h*x+c
    mu = (h*k*u/R+h*v/y)/2
    nu = (k*u+(x/2-s2)*v/(h*y))/2
    Z = s2*q*(1+x)/x.sq()
    omega = nu-Z/(4*mu)
    gain = k.sq()*u/(8*R)+2*s*v*(1+x)/(3*sqrt3*x)
    return x, y, u, v, R, mu, nu, Z, omega, gain


def omega(xi, q):
    return data(xi, q)[8]


def gain(xi, q):
    return data(xi, q)[9]


def slope(xi, q):
    """Partial derivative of omega at fixed xi (x_q=1/2)."""
    x, y, u, v, R, mu, nu, Z, om, C = data(xi, q)
    mu_q = (h*k*(-1/(2*u*R)-h*u/(2*R.sq()))
            +h*(1/(2*v*y)+v*x/(2*y*y.sq())))/2
    nu_q = (-k/(2*u)+(v/4+(x/2-s2)/(2*v))/(h*y)
            +(x/2-s2)*v*x/(2*h*y*y.sq()))/2
    Z_q = s2*((1+x)/x.sq()+q/(2*x.sq())-q*(1+x)/(x*x.sq()))
    return nu_q-Z_q/(4*mu)+Z*mu_q/(4*mu.sq())


def scalar_derivative(xi, q):
    x, y, u, v, R, mu, nu, Z, om, C = data(xi, q)
    Gamma = u*R+h*v*y
    return -(4*h+8*s)-4*Gamma*(h*u-R/u+h*y/v-h*x*v/y)


def assess(kind, box):
    if kind == 'omega':
        return omega(*box)
    if kind == 'gain':
        return gain(*box)
    if kind == 'slope':
        return slope(*box)
    if kind == 'scalar_derivative':
        return scalar_derivative(*box)
    raise AssertionError(kind)


SPECS = [
    # Refinement of the existing unique selected determinant zero.
    ('root_selector_left', 'scalar_derivative', [('239/100','239/100'), ('489/1000','489/1000')], None, '-1/100'),
    ('root_selector_right', 'scalar_derivative', [('239/100','239/100'), ('491/1000','491/1000')], '1/50', None),
    ('root_determinant_negative', 'omega', [('239/100','239/100'), ('489/1000','491/1000')], None, '-1/2000'),
    # Selector range on the two level slabs.
    ('selector_lower', 'scalar_derivative', [('21/10','21/10'), ('2/5','2/5')], None, '-1'),
    ('selector_upper_low_slab', 'scalar_derivative', [('9/4','9/4'), ('47/100','47/100')], '2/25', None),
    # Nonnegative determinant below the first working slab, and at low q.
    ('early_nonnegative', 'omega', [('21/10','21/10'), ('0','13/20')], '1/1000', None),
    ('small_population_nonnegative', 'omega', [('239/100','239/100'), ('0','1/20')], '1/100', None),
    # Upper population, until the analytic x>=4/5 screen starts.
    ('upper_population_nonnegative', 'omega', [('239/100','239/100'), ('1/2','159/200')], '1/1000', None),
    # Coarse population-gap bounds, away from the selected population.
    ('negative_part_bound', 'omega', [('239/100','239/100'), ('1/20','1/4')], '-11/200', None),
    ('gain_small', 'gain', [('21/10','239/100'), ('1/20','1/5')], '1/10', None),
    ('gain_middle', 'gain', [('21/10','239/100'), ('1/5','1/4')], '17/100', None),
    ('gain_large', 'gain', [('21/10','239/100'), ('1/4','1/2')], '19/100', None),
    # Fixed-level derivatives, not derivatives along the selected curve.
    ('slope_low_slab', 'slope', [('21/10','9/4'), ('1/4','47/100')], '0', '43/100'),
    ('slope_high_slab', 'slope', [('9/4','239/100'), ('1/4','1/2')], '0', '43/100'),
]


def certify(spec):
    name, kind, rawbox, lower, upper = spec
    initial = [(Q(a),Q(b)) for a,b in rawbox]
    low = None if lower is None else Q(lower)
    high = None if upper is None else Q(upper)
    leaves = []
    todo = [(initial,0)]
    while todo:
        box, depth = todo.pop()
        image = assess(kind, [I.box(a,b) for a,b in box])
        passes = ((low is None or Q(image.lo,S) >= low)
                  and (high is None or Q(image.hi,S) <= high))
        if passes:
            leaves.append({'box': [[str(a),str(b)] for a,b in box],
                           'image_numerators': [str(image.lo),str(image.hi)]})
            continue
        assert depth < 32, (name,box,image)
        widths = [(b-a)/(ib-ia) if ib>ia else Q(-1)
                  for (a,b),(ia,ib) in zip(box,initial)]
        axis = max(range(len(box)), key=lambda j: widths[j])
        a,b = box[axis]
        assert a < b, (name,box,image)
        mid=(a+b)/2
        left=box.copy(); right=box.copy()
        left[axis]=(a,mid); right[axis]=(mid,b)
        todo.extend([(right,depth+1),(left,depth+1)])
    result={'name':name,'expression':kind,'domain':rawbox,
            'lower':lower,'upper':upper,'leaves':leaves,
            'aggregate_numerators':[
                str(min(int(v['image_numerators'][0]) for v in leaves)),
                str(max(int(v['image_numerators'][1]) for v in leaves))]}
    print(name, 'leaves',len(leaves), 'aggregate',
          [str(Q(int(z),S)) for z in result['aggregate_numerators']], flush=True)
    return result


def verify_cover(node, leaves):
    """Independent exact volume/disjoint-interior check for rectangular leaves."""
    domain = [(Q(a),Q(b)) for a,b in node]
    active=[j for j,(a,b) in enumerate(domain) if a<b]
    boxes=[[(Q(a),Q(b)) for a,b in row['box']] for row in leaves]
    for box in boxes:
        assert all(a <= x <= y <= b for (a,b),(x,y) in zip(domain,box))
    def volume(box):
        v=Q(1)
        for j in active:
            v*=box[j][1]-box[j][0]
        return v
    assert sum(map(volume, boxes),Q(0)) == volume(domain)
    # A sweep checks disjoint interiors without assuming the generation tree.
    if not active:
        assert len(boxes)==1
        return
    axis=active[0]
    ordered=sorted(boxes,key=lambda b:b[axis][0])
    for i,box in enumerate(ordered):
        for other in ordered[i+1:]:
            if other[axis][0]>=box[axis][1]:
                break
            assert any(min(box[j][1],other[j][1]) <= max(box[j][0],other[j][0]) for j in active)


def main():
    results=[certify(spec) for spec in SPECS]
    for row in results:
        verify_cover(row['domain'],row['leaves'])
        for leaf in row['leaves']:
            im=assess(row['expression'],[I.box(a,b) for a,b in leaf['box']])
            assert [str(im.lo),str(im.hi)]==leaf['image_numerators']
    payload={'kind':'exact_outward_dyadic_continuous_certificate',
             'scale':str(S),'source_commit':'2c54932b6f9a7eb9a838e5be71e91d97e10027be',
             'local_IV_commit':'e6854b4057168ac0ffbbbc0b3e4e1468203de0f1',
             'program_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
             'certificates':results}
    dest=Path(__file__).with_name('certificate.json')
    dest.write_text(json.dumps(payload,sort_keys=True,separators=(',',':'))+'\n')
    print('certificate_sha256',sha256(dest.read_bytes()).hexdigest())
    print('PASS: all continuous leaves, threshold checks and exact cover checks')


if __name__=='__main__':
    main()

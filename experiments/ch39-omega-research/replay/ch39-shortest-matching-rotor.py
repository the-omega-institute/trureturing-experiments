"""Determine the first norm-five rotor length admitting a true-zero gauge match."""
from fractions import Fraction
from itertools import combinations, product
from pathlib import Path
from runpy import run_path
from subprocess import check_output, run
from tempfile import TemporaryDirectory

from flint import arb, ctx

root = Path(__file__).parent
reference = run_path(str(root/"ch39-grid-gauge-readout.py"))
volume = reference["volume"]
model, q, anchor, ratios = (volume[k] for k in ("gauge","q","anchor","ratios"))
words = run_path(str(root/"ch39-hurwitz-gauge-words.py"))
tolerance = arb(2)**-21
match_coord = (13174679,6162390,8401872,13953440)
match_label = 11439


def accept_witness():
    a,*imag = match_coord
    quat = (-a,*imag)
    letters = words["factor"](quat,21)
    for signs in product((-1,1),repeat=3):
        signed = (-a,*(v*t for v,t in zip(imag,signs)))
        sequence = words["factor"](signed,21)
        endpoint = (1,0,0,0)
        for j in sequence:
            endpoint = words["multiply"](endpoint,words["GENERATORS"][j])
        assert endpoint == signed
    assert Fraction(8,6*5**20) == Fraction(4,286102294921875)
    time = q("26.816214799887845")
    yukawa = model["yukawa_centers"].copy()
    yukawa[2] *= q(".9936504245091254")
    state = [(arb(i)/a)**2 for i in imag] + [y*y for y in yukawa]
    for _ in range(256):
        state = reference["advance"](state,time/256,14)
    target = model["target"](match_label,time)
    residual = [1/(r*u)-target for u,r in zip(state,ratios)]
    assert all(abs(v)<q(".00000039") for v in residual)
    print("21-letter witness","".join(map(str,letters)),"residuals",*[v.str(9) for v in residual],flush=True)


class Enclosure:
    def __init__(self):
        self.targets,_ = volume["zero_targets"](tolerance)
        xl,xh,self.logs = volume["input_bounds"]()
        self.x_mid = [(l+h)/2 for l,h in zip(xl,xh)]
        self.y_mid = [v.mid() for v in self.logs]
        domain = [l.union(h) for l,h in zip(xl,xh)]
        self.jac,rates = volume["variations"](domain,self.logs,refinements=2)
        self.rate = [v/r for v,r in zip(rates,ratios)]
        self.end = model["evolve"](self.x_mid,self.y_mid,anchor,65536,diagonal_quarks=True)[0]
        self.center_error = [sum(self.jac[i,3+j]*(v-self.y_mid[j]) for j,v in enumerate(self.logs) if j != 2)
                             /ratios[i]+arb(0,tolerance) for i in range(3)]
        self.top_sensitivity = [self.jac[i,5]/r for i,r in enumerate(ratios)]

    def rough_labels(self,coord,labels,pieces):
        a,*b = coord
        dx = [(arb(a)/v)**2-c for v,c in zip(b,self.x_mid)]
        base = [(self.end[i]+sum(self.jac[i,j]*d for j,d in enumerate(dx)))/r + self.center_error[i]
                for i,r in enumerate(ratios)]
        y = self.logs[2]-self.y_mid[2]
        lo = y.lower()
        step = (y.upper()-lo)/pieces
        kept = set()
        for j in range(pieces):
            slab = (lo+j*step).union(lo+(j+1)*step)
            values = [b+s*slab for b,s in zip(base,self.top_sensitivity)]
            time = arb(0,volume["radius"])
            for i,k in combinations(range(3),2):
                allowed = (values[k]-values[i])/(self.rate[i]-self.rate[k])
                if not time.overlaps(allowed):
                    break
                time = time.intersection(allowed)
            else:
                common = (values[0]+time*self.rate[0])*(2*time).exp()
                for n in labels:
                    target = self.targets[n][0]
                    if common.overlaps(target) and self.possible_linear(target,values,time):
                        kept.add(n)
        return kept

    def possible_linear(self,target,values,time):
        for _ in range(4):
            for i in range(3):
                allowed = (target/(values[i]+time*self.rate[i])).log()/2
                if not time.overlaps(allowed):
                    return False
                time = time.intersection(allowed)
        return True

    def central_jet(self,coord,top):
        a,*imag = coord
        y = self.y_mid.copy()
        y[2] = top
        state = [(arb(i)/a)**2 for i in imag]+[v.exp() for v in y]
        for _ in range(128):
            state = reference["advance"](state,anchor/128,14)
        tube = [v+arb(0,".0625") for v in state]
        assert all(abs(v)<q("1/2") for v in reference["taylor"](tube,1)[1])
        return reference["taylor"](state,13),reference["taylor"](tube,14)[-1]

    def possible_jet(self,n,jet,top_delta):
        coeff,rem = jet
        target = self.targets[n][0]
        time = arb(0,volume["radius"])
        for _ in range(6):
            for i,r in enumerate(ratios):
                u = rem[i]
                for term in reversed(coeff):
                    u = term[i]+time*u
                value = 1/(r*u)+self.center_error[i]+self.top_sensitivity[i]*top_delta
                allowed = (target/value).log()/2
                if not time.overlaps(allowed):
                    return False
                time = time.intersection(allowed)
        return True

    def exclude(self,coord,labels):
        # Re-center only when the top input has moved enough to benefit.
        # All finer children may reuse the same center jet and variational bounds.
        y = self.logs[2]
        stack = [(y.lower(),y.upper(),labels,0,None,None)]
        while stack:
            left,right,candidates,depth,jet,center = stack.pop()
            middle = (left+right)/2
            if depth in (0,2,5,8):
                center = middle
                jet = self.central_jet(coord,center)
            relative = left.union(right)-center
            active = [n for n in candidates if self.possible_jet(n,jet,relative)]
            if not active:
                continue
            assert depth < 13,(coord,active,left,right)
            stack.extend([(left,middle,active,depth+1,jet,center),
                          (middle,right,active,depth+1,jet,center)])


def exclude_other_orbits(enclosure):
    assert all(not words["hits"](d) for d in range(16))
    counts = (4,39,158,750,3807,19485)
    with TemporaryDirectory() as temp:
        binary = str(Path(temp)/"shell-scan")
        run(["c++","-O3","-std=c++17",str(root/"ch39-hurwitz-shell-scan.cpp"),"-o",binary],check=True)
        for d,count in enumerate(counts,16):
            rows = check_output([binary,str(d)],text=True).splitlines()
            coords = [tuple(map(int,row.split())) for row in rows]
            assert len(set(coords)) == len(coords) == count
            if d == 16:
                assert set(coords) == {tuple(map(abs,quat)) for quat,letters in words["hits"](d)}
            if d == 21:
                assert match_coord in coords
            unresolved = 0
            for coord in coords:
                labels = enclosure.rough_labels(coord,enclosure.targets,1)
                if labels:
                    labels = enclosure.rough_labels(coord,labels,32)
                if d == 21 and coord == match_coord:
                    labels.discard(match_label)
                if labels:
                    unresolved += len(labels)
                    enclosure.exclude(coord,labels)
            print("depth",d,"in-box absolute tuples",count,"fine exclusions",unresolved,flush=True)


def verify():
    ctx.prec = 128
    accept_witness()
    enclosure = Enclosure()
    exclude_other_orbits(enclosure)
    center = model["yukawa_centers"][2]
    y_ref = 2*(center*q(".9936504245")).log()
    jet = enclosure.central_jet(match_coord,y_ref)
    possible = []
    for k in range(9900,10100):
        factors = (arb(k)/10000).union(arb(k+1)/10000)
        if enclosure.possible_jet(match_label,jet,2*(center*factors).log()-y_ref):
            possible.append(k)
    assert min(possible) == 9921 and max(possible) == 9946
    fixed_time = q("26.816214799887845")
    top_window = q(".993645").union(q(".993655"))
    shift = 2*(center*top_window).log()-y_ref
    for i,r in enumerate(ratios):
        # Re-evaluate the center trajectory at the witness time; use the
        # variational bounds for all nine low-energy matter coordinates.
        coeff,rem = jet
        u = rem[i]
        for term in reversed(coeff):
            u = term[i]+(fixed_time-anchor)*u
        other_error = sum(enclosure.jac[i,3+j]*(y-enclosure.y_mid[j])
                          for j,y in enumerate(enclosure.logs) if j!=2)/r
        err = 1/(r*u)-model["target"](match_label,fixed_time)
        err += other_error + enclosure.top_sensitivity[i]*shift
        assert abs(err) < tolerance
    print("First matching length 21, eight signed words, unique attainable zero label 11439")


if __name__ == "__main__":
    verify()

"""Calibrate two angles at the least contact strength with two radial mediator copies."""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat

root = Path(__file__).parent
circle = run_path(str(root/"ch39-two-tick-cp-circle.py"))
model = run_path(str(root/"ch39-four-contact-channels.py"))
sine12, sine23 = model["s"], model["t"]
eye = arb_mat([[1,0],[0,1]])


def constraints(point, response, gain):
    f, jacobian = model["observables"]((gain,*point), True, response)
    residual = arb_mat([[f[i,0]+s*s*f[2,0]] for i,s in enumerate((sine12,sine23))])
    slopes = arb_mat([[jacobian[i,j]+s*s*jacobian[2,j] for j in (1,2)]
                      for i,s in enumerate((sine12,sine23))])
    return residual, slopes


def isolate(seed, response, gain):
    point = list(map(arb,seed))
    for _ in range(4):
        f, jac = constraints(point, response, gain)
        step = jac.solve(f)
        point = [(x-step[i,0]).mid() for i,x in enumerate(point)]
    f, jac = constraints(point,response,gain)
    inverse = jac.inv()
    offset = inverse*f
    boxes = []
    for radius in map(arb, (".0001", "1e-8")):
        box = [x+arb(0,radius) for x in point]
        remainder = eye-inverse*constraints(box,response,gain)[1]
        lipschitz = max(sum(abs(remainder[i,j]).upper() for j in range(2)) for i in range(2))
        assert max(abs(v).upper() for v in offset.entries()) + lipschitz*radius < radius
        boxes.append(box)
    return boxes


def possible_phase(strength, response, gain, angles=(sine12,sine23)):
    # W_ij = a_ij + exp(i*phi) z_ij. Eliminate cos(phi), sin(phi) from two rows.
    action = gain*response + strength*model["contact"]
    mean = action.trace().real/2
    reduced = action-mean*model["identity"]
    x,y,z = reduced[0,0].real,reduced[0,1].real,reduced[0,1].imag
    gap = (x*x+y*y+z*z).nonnegative_part().sqrt()
    plane = model["shell"]["cis"](mean)*(gap.cos()*model["identity"]+acb(0,gap.sinc())*reduced)
    projected = model["basis"]*plane*model["basis"].transpose()
    constants, coefficients = [], []
    for i,j in model["indices"]:
        a, value = model["radial"][i,j].real, projected[i,j]
        constants.append(a*a+abs(value)**2)
        coefficients.append(2*a*value)
    rows = []
    for i,s in enumerate(angles):
        v = coefficients[i]+s*s*coefficients[2]
        rows.append((v.real,-v.imag,s*s*(1-constants[2])-constants[i]))
    (a,b,c),(d,e,f) = rows
    det, cosine, sine = a*e-d*b, c*e-f*b, a*f-d*c
    if 0 not in cosine*cosine+sine*sine-det*det:
        return False, None
    if 0 in det:
        return True, None
    return True, acb(cosine/det, sine/det).arg()


def prove_least(response, gain, domain):
    pending = [(arb(0),arb("2.6"))]
    while pending:
        left, right = pending.pop()
        strength = left.union(right)
        possible, phase = possible_phase(strength,response,gain)
        if not possible:
            continue
        if phase is not None and domain[0].contains_interior(strength) and domain[1].contains_interior(phase):
            continue
        assert right-left > arb("1e-10")
        mid = ((left+right)/2).mid()
        pending.extend(((left,mid),(mid,right)))


def verify(response):
    early, _ = circle["two_ticks"]["radial_reflection"]()
    gain = 1+early
    assert arb(18)/25 in gain
    inverse = circle["verify"](response)
    circle["exclude_gain_window"](inverse, arb(0).union(arb(9)/25))
    outer, solution = isolate(("2.5605","3.0833"),response,gain)
    prove_least(response,gain,outer)
    def held_out(params):
        matrix = model["observables"]((gain,*params),response=response)
        ud,ub,cd,cb = [matrix[i,j] for i,j in ((0,0),(0,2),(1,0),(1,2))]
        product = -ud*ub.conjugate()*cd.conjugate()*cb
        invariant = product-abs(ud*cb*ub)**2/(1-abs(ub)**2)
        return abs(ub), invariant.arg()
    sine, phase = held_out(solution)
    assert arb(".003677") < sine < arb(".003787") and arb("1.141") < phase < arb("1.153")
    print("least H and common phase", *(v.str(10) for v in solution))
    print("two held-out observables",sine.str(11),phase.str(11))
    _, other = isolate(("2.6188","-2.768"),response,gain)
    second_sine, second_phase = held_out(other)
    assert second_sine > arb(".012") and second_phase < arb("1.02")
    print("higher-H branch", *(v.str(10) for v in other),second_sine.str(10),second_phase.str(10))


if __name__ == "__main__":
    verify(circle["two_ticks"]["verify"]())

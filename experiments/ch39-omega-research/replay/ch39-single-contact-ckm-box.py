"""Exclude every single incident-arc contact throughout the four CKM target intervals."""
from pathlib import Path
from runpy import run_path

from flint import arb, acb, acb_mat, arb_mat

arc = run_path(str(Path(__file__).with_name("ch39-local-flavor-arc-count.py")))
shell, fit, basis, p = (arc[k] for k in ("shell","fit","basis","p"))
cis, dot = fit["cis"], shell["dot"]
radii = [shell["sine_radius"]]*3 + [shell["phase_radius"]]
parameters = [v+arb(0,r.upper()) for v,r in zip(shell["centers"],radii)]


def data_derivatives(args, normal, derivatives=True):
    x,y,s,t,u,d = args
    a,b,c = ((1-z*z).sqrt() for z in (s,t,u))
    e = cis(d)
    g12 = acb_mat([[a,s,0],[-s,a,0],[0,0,1]])
    g13 = acb_mat([[c,0,u*e.conjugate()],[0,1,0],[-u*e,0,c]])
    g23 = acb_mat([[1,0,0],[0,b,t],[0,-t,b]])
    right = acb_mat([[cis(x),0,0],[0,cis(y),0],[0,0,1]])
    m = g23*g13*g12*right
    def left_phase(m):
        return [(sum(m[i,j]*p[j] for j in range(3))/p[i]).conjugate() for i in range(3)]
    def trace_term(l,m):
        return sum(normal[j,i]*l[i]*m[i,j] for i in range(3) for j in range(3))
    l = left_phase(m)
    f = [abs(l[0])**2-1, abs(l[1])**2-1, trace_term(l,m)]
    if not derivatives:
        return f
    dg12 = acb_mat([[-s/a,1,0],[-1,-s/a,0],[0,0,0]])
    dg23 = acb_mat([[0,0,0],[0,-t/b,1],[0,-1,-t/b]])
    dg13 = acb_mat([[-u/c,0,e.conjugate()],[0,0,0],[-e,0,-u/c]])
    i = acb(0,1)
    by_delta = acb_mat([[0,0,-i*u*e.conjugate()],[0,0,0],[-i*u*e,0,0]])
    directions = [m*acb_mat([[i,0,0],[0,0,0],[0,0,0]]),
                  m*acb_mat([[0,0,0],[0,i,0],[0,0,0]]),
                  g23*g13*dg12*right, dg23*g13*g12*right,
                  g23*dg13*g12*right, g23*by_delta*g12*right]
    cols = []
    for dm in directions:
        dl = left_phase(dm)
        cols.append((2*(dl[0]*l[0].conjugate()).real,
                     2*(dl[1]*l[1].conjugate()).real,
                     trace_term(dl,m)+trace_term(l,dm)))
    return cols


def excludes_by_derivative(args, normal):
    center = [v.mid() for v in args]
    f = data_derivatives(center,normal,False)
    slopes = data_derivatives(center,normal)
    jac = acb_mat([col[:2] for col in slopes[:2]])
    coeff = jac.solve(acb_mat([[col[-1]] for col in slopes[:2]]))
    l,m = coeff[0,0],coeff[1,0]
    value = f[2]-l*f[0]-m*f[1]
    derivatives = data_derivatives(args,normal)
    weights = []
    for x,c,(dx,dy,dq) in zip(args,center,derivatives):
        derivative = dq-l*dx-m*dy
        value += (x-c)*derivative
        weights.append(float((abs(derivative)*(x-c).rad()).upper()))
    return 0 not in value, weights


def resolve_cell(args,normal):
    pending = [(args,0)]
    count = 0
    while pending:
        box, depth = pending.pop()
        count += 1
        excluded, weights = excludes_by_derivative(box,normal)
        if excluded:
            continue
        assert depth < 22
        dim = max(range(6),key=weights.__getitem__)
        v = box[dim]
        midpoint = v.mid()
        for split in (v.lower().union(midpoint), midpoint.union(v.upper())):
            child = box.copy()
            child[dim] = split
            pending.append((child,depth+1))
    return count


def exclusion_count(index, contacts, admissible=None):
    native = arc["native"]
    t = native["plane_response"](acb.zeta_zero(index).imag)
    normals = [basis*(t*q-q*t)*basis.transpose() for q in contacts]
    return exclude_normals(normals,admissible)


def exclude_normals(normals,admissible=None):
    v = shell["target_box"]()
    offsets, coeffs = shell["rephasing_equations"](v)
    pending = [(0,0,0,normals)]
    checks = 0
    while pending:
        level,i,j,active = pending.pop()
        unit = arb.pi()/(1<<level)
        x,y = fit["around"]([(2*i+1)*unit,(2*j+1)*unit],unit)
        h,k = cis(x),cis(y)
        ds = (cis(x-y),h,k)
        if any(0 not in a+sum((b*t).real for b,t in zip(row,ds)) for a,row in zip(offsets,coeffs)):
            continue
        rows = [[z*t for z,t in zip(row,(h,k,1))] for row in v]
        l = [(dot(row,p)/p[i]).conjugate() for i,row in enumerate(rows)]
        w = acb_mat([[z*li for z in row] for row,li in zip(rows,l)])
        if admissible is not None and not admissible(w):
            continue
        keep = [n for n in active if 0 in (n*w).trace()]
        if level == 10:
            for n in keep:
                checks += resolve_cell([x,y,*parameters],n)
        elif keep:
            pending.extend((level+1,2*i+a,2*j+b,keep) for a in (0,1) for b in (0,1))
    return checks


if __name__ == "__main__":
    for index in (11431,11439,11737):
        checks = exclusion_count(index,arc["arc_projectors"]().values())
        print(index,"all ten contacts excluded; centered checks",checks)

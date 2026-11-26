"""Test the full CKM box against all spectral functions of the native flavor response."""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb

schur = run_path(str(Path(__file__).with_name("ch39-shell-ckm-schur-gain.py")))
fit = schur["fit"]
shell = fit["shell"]
geo = schur["geo"]
p = fit["p"]
n, m = schur["axes"]
dot = shell["dot"]
norm = lambda v: arb(dot(v,v)).sqrt()
second_axis = shell["cross"](p,n)
basis = acb_mat([[x/norm(n), y/norm(second_axis)] for x,y in zip(n,second_axis)])
columns = basis.transpose()*acb_mat([[a/norm(n), b/norm(m)] for a,b in zip(n,m)])
inside = schur["inside"]
origin = schur["origin"]
portals = (schur["first"], schur["second"])
source = acb_mat([[int(i==j) for j in portals] for i in inside])


def plane_response(height):
    matrix = geo["matrix"](0, geo["phases"](height))
    block = acb_mat([[matrix[i,j] for j in inside] for i in inside])
    internal = block.solve(source)
    small = acb_mat([[-matrix[origin,i]*internal[inside.index(i),k]*matrix[j,origin]
                      for k,j in enumerate(portals)] for i in portals])
    return columns*small*columns.transpose()


def enclose_native_axes():
    groups = {}
    for zero in acb.zeta_zeros(10815,1278):
        t = plane_response(zero.imag)
        diag = (t[0,0]-t[1,1]).real/2
        radius = (diag**2 + abs(t[0,1])**2).sqrt()
        assert radius > 0
        axis = (acb(diag/radius), t[0,1]/radius)
        # Bins only group intervals. Every original interval enters its group's hull.
        key = tuple(round(float(v.mid())*32) for v in (axis[0].real,axis[1].real,axis[1].imag))
        prev = groups.get(key)
        groups[key] = axis if prev is None else tuple(a.union(b) for a,b in zip(prev,axis))
    return list(groups.values())


def prove_mismatch(axes):
    v = shell["target_box"]()
    offsets, coeffs = shell["rephasing_equations"](v)
    pending = [(0,0,0,axes)]
    count = 0
    while pending:
        level,i,j,remaining = pending.pop()
        count += 1
        unit = arb.pi()/(1<<level)
        x,y = fit["around"]([(2*i+1)*unit,(2*j+1)*unit],unit)
        h,k = fit["cis"](x),fit["cis"](y)
        phases = (fit["cis"](x-y), h, k)
        if any(0 not in a + sum((b*t).real for b,t in zip(row,phases))
               for a,row in zip(offsets,coeffs)):
            continue
        rows = [[z*t for z,t in zip(row,(h,k,1))] for row in v]
        left = [(dot(row,p)/p[i]).conjugate() for i,row in enumerate(rows)]
        w = basis.transpose()*acb_mat([[l*z for z in row] for l,row in zip(left,rows)])*basis
        delta = w[1,1]-w[0,0]
        keep = []
        for d,z in remaining:
            diag = z*w[1,0]-z.conjugate()*w[0,1]
            upper = 2*d*w[0,1]+z*delta
            lower = -2*d*w[1,0]-z.conjugate()*delta
            # A 2x2 unitary has equal off-diagonal moduli in the native eigenbasis.
            if not 2*abs(diag)**2 + abs(upper)**2 + abs(lower)**2 > 8*arb(".30")**2:
                keep.append((d,z))
        if not keep:
            continue
        assert level < 16, (level,i,j)
        pending.extend((level+1,2*i+a,2*j+b,keep) for a in (0,1) for b in (0,1))
    print("1278 native eigenbases, full CKM box, off-diagonal modulus > 0.30;",count,"phase boxes")


if __name__ == "__main__":
    prove_mismatch(enclose_native_axes())

"""保留共同 724 钟时的两段相邻奇异值比全域排除。"""
from pathlib import Path
from runpy import run_path
from flint import acb, acb_mat, arb, ctx

def load(name):
    return run_path(str(Path(__file__).resolve().parent / ("ch39-" + name + ".py")))

ctx.prec = 112
graph = load("three-ray-chiral-boundary")
pump = load("native-light-projector-pump")
cert = load("moving-preparation-quark-shape")
singular_values = load("first-zero-relative-quark-spectrum")["singular_values"]
adj = pump["adjoint"]
source = cert["orthonormal_sources"](graph, adj)
frames = list(pump["rotation_bases"](graph))
inputs = [adj(frame)*source for frame in frames]
pole = [graph["soft"]["pole"]]
rays = graph["boundary_sites"]()
rows = [[i for i in range(frame.ncols()) if any(frame[j,i] != 0 for j in sites)]
        for frame, sites in zip(frames, (pole, rays, rays))]
assert all(len(row) == 1 for row in rows)
pivot = acb_mat([[m[row[0],j] for j in range(3)] for m,row in zip(inputs,rows)])

root2, soft, imag = arb(2).sqrt(), arb("6e-8"), acb(0,1)

def frob(matrix):
    return sum(abs(z)**2 for z in matrix.entries()).sqrt().upper()

def verify(charge, sign, xlo, xhi, index, limit):
    shift = [acb("1.78",".25"), acb("1.80","-.09")]
    if sign < 0:
        shift = [z.conjugate() for z in shift[::-1]]
    basis_change = pivot.inv()*acb_mat([[arb(".3"),*shift],[0,1,0],[0,0,1]])
    filter_error = arb("2e-8")*frob(basis_change)
    base, links = pump["clock_matrices"](graph, charge)
    sectors = []
    for frame, incoming in zip(frames, inputs):
        a,b = [adj(frame)*link*frame for link in links]
        n = frame.ncols()
        identity = acb_mat([[int(i == j) for j in range(n)] for i in range(n)])
        transformed = incoming*basis_change
        sectors.append((adj(frame)*base*frame, a, b, adj(a)*a+a*adj(a),
                        adj(b)*b+b*adj(b), transformed, frob(transformed), identity))
    frequency = -2*sign if charge == -4 else sign
    limit = arb(limit)
    stack = [(arb(xlo),arb(xhi),arb(0),2*arb.pi(),0)]
    accepted = 0
    while stack:
        x0,x1,y0,y1,depth = stack.pop()
        x,y = (x0+x1)/2,(y0+y1)/2
        r,s = abs(frequency)*(x1-x0)/2,(y1-y0)/2
        phase_x,phase_y = (imag*frequency*x).exp(),(imag*y).exp()
        change = root2*(r+s)+soft
        center_rows, errors = [], []
        if change < arb(".55"):
            for base,a,b,a2,b2,source_rows,weight,I in sectors:
                f,g = a*phase_x,b*phase_y
                H = base+f+adj(f)+g+adj(g)
                levels,vectors = H.eig(right=True, algorithm="rump")
                order = sorted(range(len(levels)), key=lambda k:float(levels[k].real))
                k = order[0]
                energy = levels[k].real
                gap = min((levels[j].real-energy).lower() for j in order[1:])
                assert all(z.imag.contains(0) for z in levels)
                if not gap > 2*change:
                    break
                v = acb_mat([[vectors[j,k]] for j in range(H.nrows())])
                v = v/(adj(v)*v)[0,0].real.sqrt()
                row = adj(v)
                residual = soft
                for radius,square in ((r,a2),(s,b2)):
                    residual += radius*(row*square*v)[0,0].real.max(arb(0)).sqrt()
                sine = residual/(gap-change)
                P = v*row
                reduced = (H-energy*I+P).inv()*(I-P)
                response = adj(source_rows)*reduced
                linear = sum(rad*frob(response*((z-adj(z))*v))
                             for rad,z in ((r,f),(s,g)))
                tail = sine**2*(1+residual/gap)
                tail += (2*change*sine+root2*(r*r+s*s)/2+soft)/gap
                errors.append((root2*sine*weight).min(linear+weight*tail))
                center_rows.append((row*source_rows).entries())
        if len(center_rows) == 3:
            error = sum(e**2 for e in errors).sqrt()+filter_error
            try:
                sv = singular_values(acb_mat(center_rows))
            except (AssertionError, ValueError):
                sv = None
            if sv is not None and sv[0] < sv[1] < sv[2] and sv[index] > error:
                if (sv[2]+error)/(sv[index]-error) < limit:
                    accepted += 1
                    continue
        assert depth < 24
        if r > s:
            stack.extend([(x0,x,y0,y1,depth+1),(x,x1,y0,y1,depth+1)])
        else:
            stack.extend([(x0,x1,y0,y,depth+1),(x0,x1,y,y1,depth+1)])
    print(charge, sign, xlo, xhi, index, limit, accepted, flush=True)
    return accepted

total = 0
for sign in (1,-1):
    for lo,hi,outer in ((0,".41",True),(".41","2.5",False),("2.5",arb.pi(),True)):
        tests = ((-4,1,"2.65"),(2,0,"13.9")) if outer else ((-4,0,"8.5"),(2,1,"5.45"))
        for charge,index,limit in tests:
            total += verify(charge, sign, lo, hi, index, limit)
assert arb("2.65")*arb("13.9")*(arb(9)/7)**2 < arb("61.32")
assert arb("8.5")*arb("5.45")*(arb(9)/7)**2 < arb("77.12")
print("accepted", total)
cert["target_gaps"]()

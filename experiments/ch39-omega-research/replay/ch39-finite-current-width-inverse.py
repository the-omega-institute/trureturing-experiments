"""Recover the absolute scale and site coupling from five finite ratios and a width."""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

root = Path(__file__).parent
widths = run_path(str(root/"ch39-current-width-normalization.py"))
study, four, q = widths["study"], widths["four"], widths["q"]
nodes = tuple(map(q, ("1/1024","1/50","1/5","1")))
extra_node = q("1/2")
unit, colour, right = study["ports"][-3:]


def green(e, derivative=False):
    return 3*colour(e,unit,derivative)/8 + right(e,unit,derivative)/4


def calibration_box(readings):
    neutral = widths["inverse"]["neutral"]
    _, _, cp, rp = neutral["graph_constants"](port_builder=neutral["pole_polynomials"])
    for _, denominator in (cp,rp):
        for value, multiplicity in denominator.complex_roots():
            assert multiplicity == 1 and value.imag == 0
            assert value.real == 0 or q(".04") < value.real*unit < q("3.4")
    assert four["z_limit"]+green(q(".03")) < 0 < green(q(".04"))
    xs = (arb(0),*nodes,extra_node)
    row = [arb(1),*readings]
    differences = [row[0]]
    for j in range(1,6):
        row = [(b-a)/(xs[i+j]-xs[i]) for i,(a,b) in enumerate(zip(row,row[1:]))]
        differences.append(row[0])
    mean_pole = -differences[2]/differences[3]-xs[3]
    assert mean_pole/q(".03") < q("1.4")
    x,y,z = xs[3:]
    m0 = differences[5]
    m1 = -differences[4]-z*m0
    m2 = differences[3]-(y+z)*m1-y*z*m0
    m3 = -differences[2]-(x+y+z)*m2-(x*y+x*z+y*z)*m1-x*y*z*m0
    cut = q(".02")
    assert m1-cut*m0 > 0
    assert (m1-cut*m0)*(m3-cut*m2)-(m2-cut*m1)**2 > 0
    # At least two poles exceed cut. Interlacing places the second largest
    # below the last bare graph pole, in units of the unknown scale.
    assert cut/q("3.4") > q(".0058")
    return q(".0058").union(q("1.4"))


def response(x, L, z):
    # a(-x/L) increases with L. Its L derivative decreases by positivity
    # of the Euclidean graph spectral measure, including the massless pole.
    lo, hi = L.lower(), L.upper()
    a = green(-x/lo).union(green(-x/hi))
    da = (x/hi**2*green(-x/hi,True)).union(x/lo**2*green(-x/lo,True))
    h = 1/(x*(z+a))
    return h, -x*h**2*da, -x*h**2


def equations(L, z, observed):
    columns = list(zip(*(response(x,L,z) for x in nodes),strict=True))
    rows = [[col[i]-col[3]-v*(col[0]-col[3]) for col in columns]
            for i,v in enumerate(observed,1)]
    return arb_mat([[r[0]] for r in rows]), arb_mat([r[1:] for r in rows])


def locate(readings,*,refine=False):
    v = [(r-1)/x for r,x in zip(readings[:4],nodes,strict=True)]
    observed = [(w-v[3])/(v[0]-v[3]) for w in v[1:3]]
    domain = calibration_box(readings)
    step = (domain.upper()-domain.lower())/64
    pending = [((domain.lower()+i*step).union(domain.lower()+(i+1)*step),
                arb(0).union(four["z_limit"])) for i in range(64)]
    solutions = []
    identity = arb_mat([[1,0],[0,1]])
    while pending:
        L,z = pending.pop()
        value, jac = equations(L,z,observed)
        if any(not v.contains(0) for v in value.entries()):
            continue
        centre = arb_mat([[L.mid()],[z.mid()]])
        f, j = equations(L.mid(),z.mid(),observed)
        if not j.mid().det().contains(0):
            B = j.mid().inv()
            image = centre-B*f+(identity-B*jac)*(arb_mat([[L],[z]])-centre)
            if not (L.overlaps(image[0,0]) and z.overlaps(image[1,0])):
                continue
            nl, nz = L.intersection(image[0,0]), z.intersection(image[1,0])
            if nl.rad() < q("2e-6") and nz.rad() < q("2e-6"):
                if refine and (nl.rad()<L.rad()/2 or nz.rad()<z.rad()/2):
                    pending.append((nl,nz))
                    continue
                solutions.append((nl,nz))
                continue
            if max(nl.rad()/L.rad(),nz.rad()/z.rad()) < q(".8"):
                pending.append((nl,nz))
                continue
        if L.rad() > 30*z.rad():
            pending.extend([(L.lower().union(L.mid()),z),(L.mid().union(L.upper()),z)])
        else:
            pending.extend([(L,z.lower().union(z.mid())),(L,z.mid().union(z.upper()))])
    assert solutions
    L,z = solutions[0]
    for l1,z1 in solutions[1:]:
        L,z = L.union(l1),z.union(z1)
    first, last = [response(x,L,z)[0] for x in (nodes[0],nodes[-1])]
    rho = (v[0]-v[3])/(4*L*(first-last))
    eta = v[3]-4*rho*L*last
    assert min(rho,eta)>0
    return L,z,rho,eta


def inferred_parameters():
    g, u, w, _, _, active = study["ports"][:6]
    exact = [1+(x/g**2+4/(1/(2*w**2)+active(-x)))/(18*u**2) for x in (*nodes,extra_node)]
    retuned = study["retuned"]
    _, _, (alpha,beta,family), _ = retuned["current"]["light_current"](retuned["logs"])
    width = widths["reference_width"]((alpha,beta,family))
    result = invert_observations(exact,width)
    target,L = result[:2]
    time = study["nominal"][2]
    assert L.contains(1)
    assert study["comp"]["model"]["target"](11430,time)>q("3.29412713")>target
    assert target>q("3.29402558")>study["comp"]["model"]["target"](11432,time)
    assert abs(L.log())/2 < q("3.233e-8")
    return exact,*result


def invert_observations(readings,width,relative_shape="1e-9",relative_width="1e-5",*,index=None):
    error = arb(0,q(relative_shape).upper())
    L,z,rho,eta = locate([v*(1+error) for v in readings],refine=index is not None)
    width *= 1+arb(0,q(relative_width).upper())
    c = 1/L
    f0 = widths["normalized_lepton_width"](rho,eta,c,z)/width
    offset = widths["inverse"]["port_moments"]()[-1]
    target = 15*L*(rho*f0-offset)
    if index is None:
        rho,eta,c = rho*f0,eta*f0,c*f0
        masses = four["physical_masses"](c,rho,z)
    else:
        at = study["comp"]["model"]["target"]
        time = study["nominal"][2]
        assert at(index-1,time)>target>at(index+1,time)
        masses = four["physical_masses"](c,rho,z)
        physical_rho = at(index,time)/(15*L)+offset
        f0 = f0.intersection(physical_rho/rho)
        rho,eta,c = physical_rho,eta*f0,c*f0
    return target,L,z,rho,eta,c,masses


def verify():
    exact,target,L,z,rho,eta,c,masses = inferred_parameters()
    time = study["nominal"][2]
    thresholds = run_path(str(root/"ch39-neutral-current-threshold-inverse.py"))
    thresholds["exclude_lower_tops"](*masses,c,z/rho,eta,1/rho,time+L.log()/2)
    print("five normalized inverse-current readings",*(v.str(15) for v in exact))
    print("scale log error",(abs(L.log())/2).upper().str(10))
    print("six finite momenta target bounds",target.lower().str(15),target.upper().str(15))


if __name__ == "__main__":
    verify()

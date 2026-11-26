"""Finite coherent depth for return-normalized sampling of the four-tap filter."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx

shell = run_path(str(Path(__file__).with_name('ch39-refined-mass-obstruction.py')))


def sectors(nodes, edges):
    roots = []
    visited = set()
    for index, p in enumerate(nodes):
        if index in visited:
            continue
        orbit = [index]
        for _ in range(2):
            p = p[1:]+p[:1]
            j = nodes.index(p)
            if j not in orbit:
                orbit.append(j)
        visited.update(orbit)
        roots.append(orbit)
    zeta = acb(-arb(1)/2, arb(3).sqrt()/2)
    blocks = []
    for m in range(3):
        orbits = [o for o in roots if m==0 or len(o)>1]
        columns = {j:(i,zeta**(m*r)/arb(len(o)).sqrt()) for i,o in enumerate(orbits) for r,j in enumerate(o)}
        matrices = [acb_mat(len(orbits),len(orbits)) for _ in range(3)]
        for i,j,weight,charge in edges:
            if i not in columns or j not in columns:
                continue
            a,x = columns[i]
            b,y = columns[j]
            coefficient = x.conjugate()*y*shell['ball'](weight)
            matrices[charge][a,b] += coefficient
            matrices[-charge][b,a] += coefficient.conjugate()
        blocks.append(matrices)
    return roots,blocks


def at(block,t):
    w = (1+acb(0,t))/(1-acb(0,t))
    return block[0] + block[1]*w + block[-1]*w.conjugate()


def finite_readout():
    cutoff = arb('.070')
    degree = 32
    # A positive polynomial, bounded by one after this common fixed normalization.
    normalization = (degree*(arb('.079546')/cutoff).acosh()).cosh()
    amplified = (degree*(arb('.0795438')/cutoff).acosh()).cosh()
    leakage = 1/(amplified*arb('.001'))**2
    port_floor = (amplified/normalization)**2*arb('.001')**2*(1-leakage)
    assert port_floor > arb('.000000996')
    relative_rounding = 2/(arb(2)**45*port_floor)
    retained = ((1-leakage)/(1+leakage)*(1-relative_rounding)/(1+relative_rounding))**2
    coefficients = 4*(1-retained)*(arb(13)/4)**2
    h = arb(2).sqrt()
    eta = coefficients*3*(1+h)/((3*h-4)*(9-4*h)/13)
    shift = 2*(1-arb('.99839735')-arb('8.762e-6'))*eta/(1-eta)
    lo = arb('277.636')*((1-eta)/(1+eta)).sqrt()
    hi = arb('277.641')*((1+eta)/(1-eta)).sqrt()
    assert eta < arb('.001467') and shift < arb('.00000469')
    assert 0 < arb('8.762e-6')-shift < arb('8.782e-6')+shift < 1/arb(2)**16
    assert lo > arb('277.22') and hi < arb('278.06')
    print('Degree 64 with 45-bit absolute readout:', lo.lower().str(10), hi.upper().str(10),
          'port floor', port_floor.lower().str(10), 'time half-width', shift.upper().str(10))


def verify():
    ctx.prec = 160
    nodes, edges = shell['geometry']()
    roots, blocks = sectors(nodes,edges)
    charge_rows = [Q(0)]*len(nodes)
    for i,j,c,k in edges:
        if k:
            charge_rows[i] += c
            charge_rows[j] += c
    lip = 2*shell['ball'](max(charge_rows))
    assert lip < arb('.0021')
    cells = 48
    error = lip/(56*cells)
    top_min,top_max,rest_ratio = arb('.0795'),arb('.0796'),arb('.99997')
    for j in range(cells):
        t = arb(2*j+1)/(56*cells)
        lifted = [[z.real for z in at(block,t).eig(multiple=True)] for block in blocks]
        by_size = sorted(lifted[0],key=lambda v:float(v.mid()))
        winner = by_size[-1]
        others = by_size[:-1] + lifted[1] + lifted[2]
        assert top_min < winner-error and winner+error < top_max
        assert all(abs(v)+error < rest_ratio*(winner-error) for v in others)
        assert winner-error > arb('.0795438')
        assert all(abs(v)+error < arb('.0795412') for v in others)
        assert all(abs(v)+error < arb('.070') for v in by_size[:-1])
        assert winner+error < arb('.079546')

    # Only one known eigenpair is needed for a sample-independent component bound.
    e,vectors = at(blocks[0],arb(0)).eig(right=True)
    j = max(range(len(e)),key=lambda i:float(e[i].real.mid()))
    norm = sum(abs(vectors[i,j])**2 for i in range(len(e))).sqrt()
    vector = [v/norm for v in [vectors[i,j] for i in range(len(e))]]
    lifted = {i:vector[a]/arb(len(o)).sqrt() for a,o in enumerate(roots) for i in o}
    residual = [arb(0)]*len(nodes)
    for i,j,c,k in edges:
        if k:
            c = shell['ball'](c)/14
            residual[i] += c*abs(lifted[j])
            residual[j] += c*abs(lifted[i])
    # Weyl separates the full trivial-sector complement below .070.
    sine = sum(x*x for x in residual).sqrt()/(top_min-arb('.070'))
    minimum = min(abs(lifted[i]) for i,p in enumerate(nodes[:19]) if p!=(51,51,51))
    lower = minimum-arb(2).sqrt()*sine
    assert lower > arb('.001')
    finite_readout()
    count = 1000000
    leakage = rest_ratio**count/arb('.001')**2
    assert leakage < arb('1e-7')
    # Normalized amplitudes lose at most this much in Born probability.
    defect = 1-((1-leakage)/(1+leakage))**2
    # Total variation of the triple law is at most twice the one-step bound.
    # For each k the signed coefficient observable ranges within [-M,M].
    coeff = 4*defect*(arb(13)/4)**2
    assert coeff < arb('.00002')
    scale = arb('.0795412')
    passes = 2000
    cheb_leak = 1/((passes*(arb('.0795438')/scale).acosh()).cosh()*arb('.001'))**2
    assert cheb_leak < leakage
    # Coefficient perturbations bound the modulus error of the complex filter for every zero.
    h = arb(2).sqrt()
    multiplier = (3*h-4)*(9-4*h)/13
    relative = coeff*3*(1+h)/multiplier
    assert relative < arb('.00185')
    pivot_lower = arb('277.636')*((1-relative)/(1+relative)).sqrt()
    pivot_upper = arb('277.641')*((1+relative)/(1-relative)).sqrt()
    q_upper = 1-arb('.99839735')-arb('8.762e-6')
    time_error = 2*q_upper*relative/(1-relative)
    assert pivot_lower > 277 and pivot_upper < 279
    assert 0 < arb('8.762e-6')-time_error
    assert arb('8.782e-6')+time_error < 1/arb(2)**16
    print('Chebyshev half-degree',passes,'relative filter error',relative.upper().str(10),
          'pivot range',pivot_lower.lower().str(10),pivot_upper.upper().str(10),
          'extra time error',time_error.upper().str(10))
    print('Uniform minimum top component',lower.lower().str(9),'power', count,
          'coefficient error',coeff.upper().str(10))


if __name__=='__main__':
    verify()

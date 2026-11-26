"""十维模型中全部三角共用一个刚度，真空自行选取小的剩余输运。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, arb_poly, ctx

phase = run_path(str(Path(__file__).with_name("ch39-matter-number-three-spokes.py")))
soft = phase["soft"]
pole = soft["pole"]
orbit = soft["star"]["orbit"]


def zero_current(matrix, edges):
    def polynomial(m):
        return arb_poly([v.real for v in (-m).charpoly().coeffs()])
    denominator = polynomial(matrix)
    rotation = soft["star"]["charged"]["rotation"]
    result = {}
    for i,j in edges:
        if (i,j) in result:
            continue
        moved = []
        for sign in (1,-1):
            m = acb_mat(matrix)
            m[i,j] *= acb(0,sign)
            m[j,i] = m[i,j].conjugate()
            moved.append(polynomial(m))
        # Three towers of charge -5 and one of +5. The six neutral weights
        # have zero first derivative at the flat residual connection.
        numerator = moved[0]-moved[1]
        assert all(numerator[k].contains(0) for k in (17,18,19))
        value = phase["integrate_force"](arb_poly(numerator.coeffs()[:17]),denominator)/(8*arb.pi()**2)
        for _ in range(3):
            result[i,j],result[j,i] = value,-value
            i,j = rotation[i],rotation[j]
    return arb_mat([[result[e]] for e in edges])


def setup(h):
    geo = soft["geo"]
    edges = [(i,j) for i,j,_,_ in geo["edges"] if pole not in (i,j)]
    index = {e:j for j,e in enumerate(edges)}
    outer = sorted(set(range(19))-{pole})
    faces = []
    for i,j,k in combinations(outer,3):
        if all(pair in index for pair in ((i,j),(j,k),(i,k))):
            row = [0]*len(edges)
            for pair, sign in (((i,j),1),((j,k),1),((i,k),-1)):
                row[index[pair]] = sign
            faces.append(row)
    D = arb_mat(faces)
    A = D.transpose()*D
    for i in range(len(edges)):
        A[i,i] += 1
    op = max(sum(abs(A[i,j]) for j in range(len(edges))).upper() for i in range(len(edges)))
    # For the cubic Taylor error use the absolute incidence matrix in both
    # factors, since cancellation between triangle rows cannot be assumed.
    degree_bound = 1+3*max(sum(abs(row[i]) for row in faces) for i in range(len(edges)))
    assert op <= degree_bound <= 52
    base,_,charged = soft["vector_pencil"](h)
    matrix = base+charged
    w = 1-arb(23)/25*arb(2)**-20/(h*h+soft["point"]["study"]["offset"])
    for j in orbit:
        matrix[pole,pole] += w-1
        matrix[j,j] += w-1
        matrix[pole,j] = matrix[j,pole] = -w
    g = zero_current(matrix,edges)
    v = -A.solve(g)
    norm = lambda z:(z.transpose()*z)[0,0].sqrt()
    gnorm, vnorm = norm(g), norm(v)
    # Harmonic representative of this trial phase, then a constant-vector
    # trial state in that gauge. The pole is the reference vertex.
    incidence = arb_mat([[int(j==k)-int(i==k) for k in outer] for i,j in edges])
    lap = incidence.transpose()*incidence
    for i,j in enumerate(outer):
        lap[i,i] += w if j in orbit else 1
    potential = -lap.solve(incidence.transpose()*v)
    residual = v+incidence*potential
    squared = (residual.transpose()*residual)[0,0]+sum(
        (w if j in orbit else 1)*potential[i,0]**2 for i,j in enumerate(outer))
    print("Matter force norm, linear phase norm, linear light mass",gnorm.str(10),vnorm.str(10),
          (squared/19).sqrt().str(10))
    return degree_bound, gnorm, vnorm, squared, w


def light(K, op, vnorm, squared, w):
    trial = vnorm/K
    # ||D^T(sin(Dt)-Dt)|| <= op * ||t||^3 / 2, including pole triangles.
    error = (40*trial+K*op*trial**3/2)/(arb(".999")*K-40)
    S = squared/K**2
    ten = run_path(str(Path(__file__).with_name("ch39-minimal-shifted-ten.py")))
    # After site rephasing v/K has co-closed weighted representative
    # of energy S. Temple's denominator exceeds one.
    assert min(z.lower() for z in ten["crossing"]["shell_spectra"]((0,))[0])-2*trial-S > 1
    trial_low = S/19*(1-S/(12*w))-S**2/19
    actual = (trial_low.sqrt()-error).union((S/19).sqrt()+error)
    assert actual > 0
    return actual**2,trial+error


def verify():
    ctx.prec=256
    h=acb.zeta_zero(11431).imag
    op,gnorm,vnorm,squared,w = setup(h)
    ten = run_path(str(Path(__file__).with_name("ch39-minimal-shifted-ten.py")))
    lower,upper = arb(1702000),arb(1703300)
    # Across all zero heights and phases the normalized ten-component energy
    # belongs to [-28.5,0] by the positive closed-walk bound.
    degrees = {i:0 for i in range(19)}
    for i,j,_,_ in soft["geo"]["edges"]:
        degrees[i] += 1
        degrees[j] += 1
    h1 = acb.zeta_zero(1).imag
    loss = arb(23)/25*arb(2)**-20/(h1*h1+soft["point"]["study"]["offset"])
    dlow = {i:d-loss*(3 if i==pole else int(i in orbit)) for i,d in degrees.items()}
    bound = 10*4*18*sum(1/(dlow[i]*dlow[j]).sqrt() for i,j,_,_ in soft["geo"]["edges"])/(24*arb.pi()**2)
    assert bound < arb("28.5")
    assert arb("5.2e-19")*bound < arb("1.5e-17")
    radius = arb.pi()*(bound/(2*lower)).sqrt()
    assert (arb(3).sqrt()*radius).cos() > arb(".999")
    # ||f'||<=300 and the positive divided differences of f' are <=11
    # for f(m)=m^2 log m. Negative contributions cost <=8/e per block.
    assert 10*(600+22)/(16*arb.pi()**2) < 40
    at_low, norm_bound = light(lower,op,vnorm,squared,w)
    at_high, _ = light(upper,op,vnorm,squared,w)
    box = arb("5.8317e-19").union(arb("5.8319e-19"))
    assert at_low > box and at_high < box
    # The unique global vacuum is continuous in K. Select the segment
    # connecting the two soft-mass faces on which the mass stays in the box.
    mu = arb("6.3573e15").union(arb("6.3575e15"))
    common, mass = ten["matching_box"](box,arb(0,(arb(2).sqrt()*norm_bound).upper()),mu)
    assert arb("53.9033") < common < arb("53.9035")
    assert arb("5.687e5") < mass < arb("5.689e5")
    print("Stiffness bracket, light eigenvalue on its faces",
          lower.str(8),upper.str(8),at_low.str(9),at_high.str(9))
    print("Global ten vacuum inverse, light mass GeV, feedback guarantee",
          common.str(10),mass.str(10),arb("5.2e-19").str(4))


if __name__ == "__main__":
    verify()

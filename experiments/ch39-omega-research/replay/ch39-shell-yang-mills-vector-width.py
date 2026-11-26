"""等耦合站点 Yang–Mills 模型中，中性重矢量向最低带荷矢量对的部分宽度。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

neutral = run_path(str(Path(__file__).with_name("ch39-neutral-vector-three-body-cut.py")))
soft = neutral["shell"]["soft"]


def polarization_sum(x):
    # Use the parent mass as unit and all momenta outgoing at the vertex.
    mass = x.sqrt()
    momentum = (arb(1)/4-x).sqrt()
    energy = arb(1)/2
    momenta = ((-1,0,0,0), (energy,0,0,momentum), (energy,0,0,-momentum))
    ex, ey = (0,1,0,0), (0,0,1,0)
    parent = (ex, ey, (0,0,0,1))
    first = (ex, ey, (momentum/mass,0,0,energy/mass))
    second = (ex, ey, (momentum/mass,0,0,-energy/mass))
    def dot(a,b):
        return a[0]*b[0]-sum(a[i]*b[i] for i in (1,2,3))
    total = arb(0)
    for eps in product(parent, first, second):
        amplitude = arb(0)
        for i,j,k in ((0,1,2),(1,2,0),(2,0,1)):
            amplitude += dot(eps[i],eps[j])*(dot(momenta[i],eps[k])-dot(momenta[j],eps[k]))
        total += amplitude**2
    return total/3


def verify():
    ctx.prec = 256
    h = acb.zeta_zero(11431).imag
    base, unbroken, broken = soft["vector_pencil"](h)
    r = soft["geo"]["theta"](804)/soft["geo"]["theta"](12)
    near_mass = neutral["invariant_masses"]()[0]
    # The first invariant positive eigenmode of the unit graph is explicit.
    pole = soft["pole"]
    unit = [int(i in (1,2,5))-int(i in (0,3,4)) for i in range(19)]
    for i in range(19):
        assert sum((base[i,j]+unbroken[i,j])*unit[j] for j in range(19)) == 11*unit[i]
    decrement = arb(23)/25*arb(2)**-20/(h*h+soft["point"]["study"]["offset"])
    for j in soft["star"]["orbit"]:
        base[pole,pole] -= decrement
        base[j,j] -= decrement
        base[pole,j] += decrement
        base[j,pole] += decrement
    eigen, vec = (base+broken).eig(right=True)
    idx = min(range(19), key=lambda j: float(eigen[j].real))
    lowest = eigen[idx].real
    assert 0 < lowest < 2 and all(eigen[j].real > 7 for j in range(19) if j != idx)
    density = [abs(vec[i,idx])**2 for i in range(19)]
    total = sum(density)
    density = [v/total for v in density]
    excess = 19*sum(v*v for v in density)-1
    assert arb(".04101967") < excess < arb(".04101970")
    # Project the perturbed invariant eigenvector onto the orthogonal
    # complement of the unperturbed simple 11 eigenspace (gap > 0.75).
    err = 4*decrement
    assert err < arb("1e-12")
    norm_error = arb(2).sqrt()*err/(arb(".75")-err)
    overlap = sum(v*n for v,n in zip(density,unit))/arb(6).sqrt()
    coupling = arb(19).sqrt()*abs(overlap+arb(0,norm_error.upper()))
    assert arb(".04788320") < coupling < arb(".04788323")
    x = lowest/(r*near_mass)**2
    beta = (1-4*x).sqrt()
    assert arb(".54629") < beta < arb(".54631")
    spins = polarization_sum(x)
    assert (spins-beta**2*(1+20*x+12*x*x)/(12*x*x)).contains(0)
    # Six complex roots with |hypercharge| = 5/sqrt(60).
    width_ratio = arb(5)/2 * coupling**2 * beta/(16*arb.pi()) * spins
    assert arb(".00024562") < width_ratio < arb(".00024565")
    daughter = lowest.sqrt()/r
    threshold = near_mass+2*daughter
    assert arb(".71406184") < threshold < arb(".71406187")
    # B and the neutral vector polarized along x, both charged daughters
    # along y: the contact vertex has magnitude 2 in group/coupling units.
    phase_space = (daughter**2*near_mass/threshold**3).sqrt()/(64*arb.pi()**2)
    spectral_lower = 2*coupling**2*(arb(5)/12)*phase_space/threshold**2
    assert spectral_lower > arb("9.9738e-7")
    indices = (arb(5)/2,arb(3)/2,arb(1))
    p = (115*indices[0]-333*indices[1]+218*indices[2])/218
    j = 5*indices[0]-3*indices[1]-2*indices[2]
    assert (p-arb(3)/109).contains(0) and j == 6
    assert p*arb(".4304") < arb(".03318") and j*arb(".4304") < arb("2.583")
    print("Relative site coupling and total squared heavy overlaps", coupling.str(12),excess.str(12))
    print("Parent and daughter masses / matching scale", near_mass.str(12),(lowest.sqrt()/r).str(12))
    print("Six-channel partial width / (diagonal coupling squared times parent mass)",width_ratio.str(12))
    print("Charged-pair-plus-neutral threshold and one-spin spectral lower coefficient",
          threshold.str(12),spectral_lower.str(12))
    return daughter, coupling, density, base+unbroken


if __name__ == "__main__":
    verify()

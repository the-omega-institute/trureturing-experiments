"""将七三角端口的正反取向同时纳入能谱，用原端口三角形比较两者。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

point = run_path(str(Path(__file__).with_name("ch39-port-coupling-competition.py")))
torus = point["torus"]
study = point["study"]
adjoint = point["adjoint"]
zmax = arb(".007335")
max_gain = arb(5)/16
floor_D = arb("1.3")-4*arb.pi()/32


def endpoints(height,vector):
    D,J = point["internal"](height)
    b = J*vector
    X = D.solve(b)
    # A concave score -b* (D+zI)^(-1) b lies above its chord. Its maximum
    # excess is bounded by zmax^2 times the second derivative bound / 8.
    error = zmax**2*sum(abs(q)**2 for q in X.entries())/(4*floor_D)
    vals = [-(adjoint(b)*point["solve"](D,b,z))[0,0].real for z in (arb(0),zmax)]
    return vals,error,sum(abs(q)**2 for q in X.entries())


def verify():
    ctx.prec = 192
    center = acb.zeta_zero(11431).imag
    vec = point["prepare"]()
    other = vec.conjugate()
    ref,center_excess,_ = endpoints(center,vec)
    c = study["offset"]
    inverse_ref = 1/(center*center+c)
    # Original marked point, nearer 724 endpoint, farther 804 endpoint.
    i = point["native"]["origin"]
    j,k = point["native"]["portals"]
    unit = point["geo"]["unit"]
    reference_triangle = unit(i,j)*unit(j,k)*unit(i,k).conjugate()
    rho = reference_triangle.imag*(vec[1,0]/abs(vec[1,0])).imag
    assert rho > arb(".098")

    torus_budget = arb(".000210")
    for z,target in zip((arb(0),zmax),ref):
        nodes,rows = torus["response_grid"](z)
        d0,_ = rows[0][0]
        vals = [[d/d0*((adjoint(other)*T*other)[0,0].real-target+torus_budget)
                 for d,T in row] for row in rows]
        torus["certify"](torus["torus_coefficients"](vals),stationary_origin=False)

    def ends(h):
        vals,excess,norm = endpoints(h,other)
        return [a-b+inverse_ref-1/(h*h+c) for a,b in zip(vals,ref)],excess,norm

    low_bound = arb(1)/100
    for zeta in acb.zeta_zeros(1,700):
        vals,_,_ = ends(zeta.imag)
        low_bound = min(low_bound,*(v.lower()-center_excess for v in vals))
    high = acb.zeta_zero(701).imag
    tail_bound = -torus_budget-center_excess+inverse_ref-1/(high*high+c)
    lower = min(low_bound,tail_bound)
    assert lower > -rho/448
    self_values,_,_ = ends(center)
    assert min(q.lower() for q in self_values)-center_excess > arb(".025")

    witness = acb.zeta_zero(88083).imag
    vals,excess,norm = ends(witness)
    ceiling = max(v.upper() for v in vals)+excess
    # If the negative-orientation Schur value at e* is negative, the true
    # energy advantage is at least its magnitude / the derivative upper bound.
    advantage = -ceiling/(1+max_gain**2*norm)
    assert advantage > rho/512
    print("Comparator splitting per unit coefficient",rho.str(14),
          "negative-sector advantage upper in units of g^2",(-lower).str(12))
    print("Zero 88083 advantage lower in units of g^2",advantage.str(12),
          "uniform normalized penalty thresholds", (advantage/rho).str(12),
          (-lower/rho).str(12))


if __name__=="__main__":
    verify()

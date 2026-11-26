"""Local gauge-input volume densities for an exact zero versus a continuous height."""
from pathlib import Path
from runpy import run_path
from flint import acb, arb, ctx

model = run_path(str(Path(__file__).with_name("ch39-floating-scale-winner.py")))
q = model["q"]


def verify():
    ctx.prec = 144
    loop, ratios = model["loop"],model["ratios"]
    beta,b,d = [model["running"][k] for k in ("BETA","B","D")]
    logs = list(map(q,("-25.4446915838 -12.9928374813 -1.5528331663 -23.8189018241 "
                       "-17.8423094899 -10.1208976879 -25.6505068564 -14.9499484326 "
                       "-9.2850641599").split()))
    for n,time,half_width in [(10871,q("26.859794675"),q(".001")),
                              (10854,q("26.861182485"),q(".001011"))]:
        jacobian_log = arb(0)
        def divergence(step, gauge, traces):
            nonlocal jacobian_log
            gauge_part = 2*sum(b[i][i]*value**2 for i,value in enumerate(gauge))/loop**2
            yukawa_part = (9*(traces[0]+traces[1])+5*traces[2])/loop
            jacobian_log += step*(gauge_part+yukawa_part)
        common = model["target"](n,time)
        inverse = [common*r for r in ratios]
        start,y,_ = model["evolve"](inverse,logs,-time,131072,diagonal_quarks=True,observe=divergence)
        assert all(abs(1/(v.sqrt()*c)-1)<half_width for v,c in zip(start,model["centers"]))
        assert all(q(".99")<(v/2).exp()/c<q("1.01") for v,c in zip(y,model["yukawa_centers"]))
        traces = [sum(v.exp() for v in logs[j:j+3]) for j in (0,3,6)]
        velocities = []
        for i,r in enumerate(ratios):
            correction = sum(bij/v-dij*t for bij,dij,v,t in zip(b[i],d[i],inverse,traces))
            velocities.append(2*common-2*(beta[i]+correction/loop)/(r*loop))
        assert 0<velocities[0]<velocities[1]<velocities[2]
        determinant = jacobian_log.exp()
        factor = determinant*ratios[0]*ratios[1]*ratios[2]
        discrete = 4*factor*sum(velocities)
        height = acb.zeta_zero(n).imag
        offset = q(model["fixed"]["OFFSET"])
        height_rate = 2*height*common/(height**2+offset)
        continuous = 4*factor*height_rate*(velocities[2]-velocities[0])
        ratio = height_rate*(velocities[2]-velocities[0])/sum(velocities)
        assert q("4e-6")<ratio<q("5e-6")
        assert (continuous/discrete).overlaps(ratio)
        print(n,"Jacobian",determinant.str(9),"discrete and continuous tube coefficients",
              discrete.str(9),continuous.str(9),"density ratio",ratio.str(11))


if __name__=="__main__":
    verify()

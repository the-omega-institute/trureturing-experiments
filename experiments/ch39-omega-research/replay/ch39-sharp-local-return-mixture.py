"""用三步返回与和乐相位的有向分量直接判别无中心零点的混合。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, arb_poly

sharp = run_path(str(Path(__file__).with_name("ch39-sharp-gauge-flavor-variance.py")))
study = sharp["study"]
g2 = acb.zeta_zero(11431).imag**2
a = study["offset"]/g2
k = (study["offset"]-arb(1)/4)/(g2+arb(1)/4)
slope = (1-a)/2+k


def return_cost(x):
    return x*x/(1-k*x)


def directed_ratio(x):
    return ((1+x)*(1-a*x)).sqrt()/(1-k*x)


def directed_cost(x):
    return 8*(1+slope*x-directed_ratio(x))


def directed_derivative(x):
    square = (1+x)*(1-a*x)
    deriv = 1-a-2*a*x
    return 8*(slope-deriv/(2*square.sqrt()*(1-k*x))-k*square.sqrt()/(1-k*x)**2)


def verify_convexity():
    # 变量 X=1+x 自零到第一枚零点的值。二阶导数同号的多项式在整段严格为负。
    high = study["center"]/(acb.zeta_zero(1).imag**2+study["offset"])
    x = arb_poly([-1,high])
    f = (1+x)*(1-a*x)
    d = 1-k*x
    curvature = -(1+a)**2*d**2+4*k*f*(1-a-2*a*x)*d+8*k*k*f**2
    from math import comb
    degree = curvature.degree()
    assert all(sum(curvature[j]*comb(i,j)/comb(degree,j) for j in range(i+1))<0
               for i in range(degree+1))
    assert 1-k*(high-1)>0 and 1-a*(high-1)>0


def verify_quarter_turn(weights):
    from fractions import Fraction
    local = run_path(str(Path(__file__).with_name("ch39-three-step-zero-readout.py")))
    nodes, edges = local["shell"]["geometry"]()
    neighbors = [{} for _ in nodes]
    for i,j,weight,charge in edges:
        neighbors[i][j] = weight, charge
        neighbors[j][i] = weight, -charge
    port = nodes.index((33,33,75))
    constant, gain, counts = local["triangle_terms"](neighbors,port)
    assert counts == {0:12,1:4} and gain == Fraction(3517425,42342346661888)
    def observed(phase):
        matrix = local["shell"]["matrices"](nodes,edges,phase)
        reading = (matrix**3)[port,port]
        assert reading.imag.contains(0)
        return reading.real
    plus,minus = observed(acb(1)),observed(acb(-1))
    assert (plus-minus).overlaps(2*local["ball"](gain))
    gaps = []
    for n in (11431,*sharp["support"]):
        zero = acb.zeta_zero(n)
        w = (zero-1)/zero
        contrast = observed(-acb(0,1)*w)-observed(acb(0,1)*w)
        y = contrast/(plus-minus)
        assert abs(y-4*zero.imag/(4*zero.imag**2+1)) < arb("1e-32")
        gaps.append(contrast)
    minimum = gaps[0]-sum(p*r for p,r in zip(weights,gaps[1:]))
    print("Sharp gap in two raw quarter-turn returns",minimum.str(16),flush=True)
    assert minimum > arb("1.27e-14")


def return_at_finite_gauge_precision():
    labels = sharp["support"]
    heights = [acb.zeta_zero(n).imag for n in labels]
    gaps = [study["center"]/(h*h+study["offset"])-1 for h in heights]
    # 精确消去普通三步返回的残差，另外两个有理权重取严格在味窗内部的构造。
    coordinate = [x/(1-k*x) for x in gaps]
    middle, third = arb(".30627112"),arb(".05349231")
    first = -(middle*coordinate[1]+third*coordinate[2]+(1-middle-third)*coordinate[3])/(coordinate[0]-coordinate[3])
    weights = [first,middle,third,1-middle-third-first]
    assert all(w>0 for w in weights)
    mean = sum(w*x for w,x in zip(weights,gaps))
    assert -arb("1.7e-9") < mean < 0
    observations = [sharp["row"](acb.zeta_zero(n))[1] for n in labels]
    A,B,U,V,J = [sum(w*y for w,y in zip(weights,col)) for col in zip(*observations)]
    real = ((1-U)*(U-V)-U*(A+B)+(1+U)/(1-U)*A*B)/2
    values = [(A/(1-U)).sqrt(),(B/(1-U)).sqrt(),U.sqrt(),acb(real,J).arg()]
    shell = study["flavor"]["shell"]
    errors = [abs(v-c)/r for v,c,r in zip(values,shell["centers"],[shell["sine_radius"]]*3+[shell["phase_radius"]])]
    assert all(e<1 for e in errors)
    target = run_path(str(Path(__file__).with_name("ch39-floating-scale-winner.py")))["target"]
    central = target(11431,arb("26.81676714125598"))
    assert 2*(arb(2)**-21)/central < arb("3e-7")
    print("Exact ordinary return, gauge mean relative shift",mean.str(15),
          "absolute shift", (mean*central).str(13),flush=True)


def verify():
    return_at_finite_gauge_precision()
    verify_convexity()
    rows = [sharp["row"](acb.zeta_zero(n)) for n in sharp["support"]]
    weights,reference = sharp["attaining_distribution"](rows,[],[])
    verify_quarter_turn(weights)
    duals = []
    for name,cost in (("return",return_cost),("directed",directed_cost)):
        objective = [cost(x) for x,_ in rows]
        minimum = sum(w*y for w,y in zip(weights,objective))
        relative = k*(1+k)*minimum if name=="return" else minimum/8
        pure = 1/(4*g2+1) if name=="return" else 4*g2.sqrt()/(4*g2+1)
        print(name,"relative separation",relative.str(18),"absolute", (pure*relative).str(18),flush=True)
        dets = []
        for U in (sharp["lo"],sharp["hi"]):
            tests = sharp["phase_tests"](U,reference)
            matrix = arb_mat([[1,x,*tests(o)] for x,o in rows])
            dets.append(matrix.det())
            c = matrix.solve(arb_mat([[z] for z in objective]))
            assert c[2,0]>0 and c[3,0]>0
            def linear(observations,c=c,tests=tests):
                first,phase = tests(observations)
                return c[2,0]*first+c[3,0]*phase
            if name=="directed":
                loss = minimum/8-abs(slope-c[1,0]/8)*arb("3e-7")
                assert loss>arb(".00084825")
            upper = sharp["compare_torus"](c[3,0],linear)
            for n in (7999,17001):
                h = acb.zeta_zero(n).imag
                x = study["center"]/(h*h+study["offset"])-1
                derivative = x*(2-k*x)/(1-k*x)**2 if name=="return" else directed_derivative(x)
                assert (derivative<c[1,0] if n==17001 else derivative>c[1,0])
                assert cost(x)-c[1,0]*x > minimum+upper,(name,upper)
            duals.append((name,cost,minimum,c[1,0],linear))
        assert dets[0]*dets[1]>0
    least = {}
    for start in range(8000,17001,256):
        count = min(256,17001-start)
        zeros = acb.zeta_zeros(start,count)
        assert len(zeros)==count
        for n,z in enumerate(zeros,start):
            if n==11431 or n in sharp["support"]:
                continue
            x,values = sharp["row"](z)
            for name,cost,minimum,b,linear in duals:
                slack = (cost(x)-minimum-b*x-linear(values)).lower()
                assert slack>0
                if name not in least or slack<least[name]:
                    least[name]=slack
    for name,slack in least.items():
        print(name,"finite dual slack",slack.str(10))


if __name__ == "__main__":
    study["torus_bound"]()
    verify()

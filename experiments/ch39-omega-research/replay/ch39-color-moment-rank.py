"""Minimum hidden-state count for four circle harmonics with tiny symmetric steps."""
from pathlib import Path
from runpy import run_path
from fractions import Fraction as F
from flint import arb, acb, arb_mat, ctx

ROOT = Path(__file__).resolve().parent
two = run_path(str(ROOT/"ch39-coarse-color-harmonics.py"))
q = two["q"]
POWERS = [0, 1, 4, 9, 16]
EPS = F(1, 2**26)
N = 2**56
# Rational nodes for sign tests (10^-18) and centers of the four existence intervals (10^-51).
SEEDS = {
    3200: (996824712649384593, 997098992373363857, 996817874100489621625830945496419673318873949507075, 997093389544986291299714858303814514883357828182464, 417123293731035275182708899689600567763300628303869, 582868507058388409393639954720688655825782806769124),
    11485: (706181769591215119, 881761873634856452, 664044801812562275782429829695444047391727044977219, 874515680094434607966190554391488818711610185039253, 73096249365888992856541107687755875919718016355192, 912567634773236597606868209962075262345286054137056),
    11486: (706127864668596278, 881747702230578841, 663989993074106455074381736556104348861441880476809, 874500885464700552273905718624164119629715070105002, 73107501876815793059894580765168471006383206570492, 912556771175569683295945226248355239013295097932027),
    11487: (706045570801900116, 881726067129517750, 663906325715118136667731038085279557740571386128882, 874478299194320161975899365283795802692925292794975, 73124692422678856456439998318528478179273355376873, 912540176040899907556827409526300710149185623557452),
    11488: (706010160542288121, 881716757623396887, 663870326802980330433323639410903225381133193508982, 874468580512108983273088988160538292288120650198552, 73132093780789992885505843491034563010379482531223, 912533031506000927890783176151378662654093697874515),
    11489: (705935362563697384, 881697092661852193, 663794290040688139985449358054139514088436125971406, 874448051441149977210345611872232964502230212390191, 73147736619337912158331098800278550693299447477979, 912517932394174057166035833283066105177969338712267),
    11490: (705900026382581418, 881687802386711925, 663758370822037062232774015799520367059323728477692, 874438353038099042565665585362231296150001208482926, 73155130771828647780338560399667993464115260926484, 912510795691661612443754503048963160347298995429740),
}


def row(value, derivative):
    if derivative:
        return [k*value**(k-1) if k else q(0) for k in POWERS]
    return [value**k for k in POWERS]


def require_three(target, nodes):
    a, b = nodes
    h = EPS**2/2
    errors = [q(F(k*(k-1), 2)*h/(1-k*h)) for k in POWERS]
    for constraints in [[(a,0),(b,0),(b,1),(q(1),0)],
                        [(q(0),0),(a,0),(b,0),(q(1),0)],
                        [(q(0),0),(a,0),(a,1),(b,0)]]:
        rows = [row(value, derivative) for value, derivative in constraints]
        coeff = arb_mat([r[:-1] for r in rows]).solve(arb_mat([[-r[-1]] for r in rows]))
        poly = [*coeff.entries(), q(1)]
        margin = sum(c*m-abs(c)*e for c,m,e in zip(poly,[q(1),*target],errors))
        assert margin > 0


def three_states(target, center):
    # Value, first derivative for four harmonics of n independent lazy steps.
    x = q(EPS/3)
    phase = acb((1-x*x)/(1+x*x),2*x/(1+x*x))
    ratios = [(1-(phase**j).real)/(1-phase.real) for j in range(1,5)]
    def curve(v, r):
        root = v**(q(1)/N)
        base = 1-r*(1-root)
        val = base**N
        return val, r*root*val/(v*base)
    def evaluate(point):
        a,b,p,s = point
        residual, jac = [], []
        for r,m in zip(ratios,target):
            fa, da = curve(a,r)
            fb, db = curve(b,r)
            residual.append([1+p*(fa-1)+s*(fb-1)-m])
            jac.append([p*da,s*db,fa-1,fb-1])
        return arb_mat(residual), arb_mat(jac)
    radius = q("1e-33")
    box = [u+arb(0,radius) for u in center]
    assert q("1/2") < box[0] < box[1] < 1 and box[2] > 0 and box[3] > 0 and box[2]+box[3] < 1
    value, jac0 = evaluate(center)
    _, jac = evaluate(box)
    inverse = arb_mat(4,4,[v.mid() for v in jac0.inv().entries()])
    assert not inverse.det().contains(0)
    diff = arb_mat([[int(i==j) for j in range(4)] for i in range(4)])-inverse*jac
    defect = inverse*value
    for i in range(4):
        assert abs(defect[i,0])+radius*sum(abs(diff[i,j]) for j in range(4)) < radius


def verify():
    ctx.prec = 450
    mu = 1000*q(F(429,16)).exp()
    power = (two["geo"]["BOOK_PLANCK"]/mu)**2/(8*arb.pi())
    for index, seeds in SEEDS.items():
        n, degree = (q(144196081),1650) if index==3200 else (power,190)
        target = two["harmonics"](acb.zeta_zero(index).imag,n,degree,4)
        nodes = [q(F(x,10**18)) for x in seeds[:2]]
        require_three(target,nodes)
        center = [q(F(x,10**51)) for x in seeds[2:]]
        three_states(target,center)
        print("Zero, disjoint-band cut points, third state weight:",
              index, nodes[0].str(9), nodes[1].str(9), (1-center[2]-center[3]).str(9))


if __name__ == "__main__":
    verify()

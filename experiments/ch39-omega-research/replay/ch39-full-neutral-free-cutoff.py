"""固定同一零点、能标和门户，用共有的有符号截止求破缺背景的完整中性匹配。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx


def response_rows(data, masses, signs, rho, portal):
    """三项 beta 依次为 QQ、QT3、色。"""
    top = arb(".8616")**2
    weak = arb(".63811")**2/2
    neutral = weak + 3*arb(".467453")**2/10
    fermion, scalar, _, gauge = range(4)
    _, _, light = data["species"][0]
    em = lambda beta: 5*beta[0]/3 + beta[1]
    top_em, top_mix, top_colour = arb(16)/9, arb(2)/3, arb(2)/3
    rows = [
        (fermion, arb(0), (em(light)-top_em, light[1]-top_mix, light[2]-top_colour)),
        (fermion, top*rho, (top_em, top_mix, top_colour)),
        (gauge, weak*rho, (-arb(22)/3, -arb(22)/3, 0)),
        (scalar, weak*rho, (arb(1)/3, arb(1)/6, 0)),
        (gauge, arb(0), (0, 0, -11)),
    ]
    layers = iter(zip(masses, signs, strict=True))
    for kind, mass, beta in data["species"][3:]:
        if kind == scalar:
            expected, sign = next(layers)
            assert expected.overlaps(mass)
            mass += sign*portal*rho
        rows.append((kind, mass, (em(beta), beta[1], beta[2])))
    assert next(layers, None) is None
    return rows, ((-12,top), (6,weak), (3,neutral))


def weak_column(x, p0, p1):
    # 每单位截止的弱差通过两个解析包围接续。小截止端由原始 E1 核控制。
    if x < p0:
        radius = 2*x*(1/p0+1/p1)
        return arb(0,radius.upper())
    reference = arb("-.074256577194663365546").union(arb("-.074256577194663365544"))
    return reference + arb(0,(arb("2.398e-21")+arb("7.468e-36")/x).upper())


def verify(root):
    gravity = run_path(str(root/"ch39-three-index-free-atoms.py"))
    column, _ = gravity["spectra"](root)
    gauge = run_path(str(root/"ch39-all-field-signed-crossing.py"))
    rotation = run_path(str(root/"ch39-rotation-portal-top-energy.py"))
    signs, masses = rotation["rotation_signs"](root, 2_306_888_100_085)
    scale = "1.176397e14"
    data = gauge["joint_system"](root, scale, order=12, last=36)
    ctx.prec = 320
    rho = (arb("246.22")/arb("1.22089e19"))**2/2
    portal = arb(".42995401872")
    p0,p1 = data["p0"],data["p1"]
    species, light = response_rows(data,masses,signs,rho,portal)
    heat = data["heat"]
    rule = gauge["angular_rule"](18,60)
    logratio = (p1/p0).log()
    pi = arb.pi()
    def responses(x):
        gauge_terms = [arb(0),arb(0),arb(0)]
        for kind, mass, beta in species:
            term = gauge["angular_exposure"](kind,mass,x,p0,p1,logratio,rule)/(4*pi)
            for i,b in enumerate(beta):
                gauge_terms[i] -= b*term
        slope = sum(g*k*x*heat(2,k*rho/x)/4 for g,k in light)
        curvature = -sum(g*k*k*heat(1,k*rho/x)/4 for g,k in light)
        for mass, sign in zip(masses,signs,strict=True):
            shifted = mass+sign*portal*rho
            slope += sign*portal*x*heat(2,shifted/x)
            curvature -= portal**2*heat(1,shifted/x)
        photon, mixed, colour = gauge_terms
        return [slope/(8*pi**2),curvature/(8*pi**2),photon-8*colour/3,mixed-colour,weak_column(x,p0,p1)],colour

    initial1,initial2,initial3 = data["initial"]
    photon0 = 5*initial1/3+initial2
    physical = [0,arb(".258"),8*initial3/3-photon0,initial3-initial2,0]
    sites = list(map(arb,("1e-50","4e-15","1e-12","3.5e-11","1e-7","1.23e-6",".5","1",
                         "1.7033e-9","2.54e-7","8.96948e-6",".0303688")))
    columns, colours = [],[]
    for x in sites:
        rows,colour = responses(x)
        columns.append(column(x) + rows)
        colours.append(colour)
    def match(indices,selected,target):
        mat = arb_mat([[columns[i][j] for j in selected] for i in indices]).transpose()
        weights = gravity["solve"](mat,arb_mat([[v] for v in target]))
        colour = initial3+sum(weights[k]*colours[i] for k,i in enumerate(indices))
        negative = sum(-w for w in weights if w<0)
        assert all(not w.contains(0) for w in weights)
        print("位置数、负量、色逆耦合、各权",len(indices),negative.str(12),colour.str(12),
              [w.str(10) for w in weights],flush=True)
        return weights,negative,colour
    weights8, cost8, colour8 = match(range(8),[0,1,2,7,8,9,10,11],
                                     [1,1,arb("1e-122"),*physical])
    assert arb("23.22418") < cost8 < arb("23.22419") and arb("38.06784") < colour8 < arb("38.06785")
    weights12, cost12, colour12 = match(range(12),range(12),
                                        [1,1,arb("1e-122"),0,0,0,0,*physical])
    assert arb("2.88377e7") < cost12 < arb("2.88378e7") and 75000 < colour12 < 75400
    return sites, weights8, weights12, species, data


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)

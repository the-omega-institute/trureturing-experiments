"""两项算术读数的预算放宽多少才可避开 11431。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb

control = run_path(str(Path(__file__).with_name("ch39-arbitrary-zero-control.py")))
free = control["free"]
bands = (
    # 积分范数界，规范与端口乘子，以及核中附加的范数系数
    ("0", ".0032", "-1.759", ".463", "0"),
    (".0032", ".0045", "-1.790", ".470", "35"),
    (".0045", ".0065", "-1.634", "1.843", "1496"),
    (".0065", ".010", "-.692", "3.497", "2612"),
    (".010", ".016", "1.308", "5.896", "3782"),
)


def exclude(data):
    values, matrices = data[-2:]
    axes = (1,5,29)
    lens = [((1-arb(axes[i])**2/867)*(1-arb(axes[j])**2/867)).sqrt()/abs(v)
            for (i,j),v in zip(control["entries"],values)]
    tests = control["base_tests"](matrices)
    scores, bonuses = [], {}
    weights = list(map(arb, (".285", ".115", ".600")))
    ua, up = arb(".000055"), arb(".006")
    phase = sum((-acb(0,s)*m for s,m in zip((1,1,-1,-1),matrices)),acb_mat(2,2))
    for band, (low,high,g,h,bonus) in enumerate(bands):
        L = arb(high)
        tilts = [min(v.upper(),(L*lens[i]).upper())/(2*(1-L*lens[i]))
                 for i,v in enumerate(map(arb,(".04",".006")))]
        for signs in product((-1,1),repeat=2):
            name = f"{band} {signs}"
            matrix = weights[2]*phase/up
            for i,sign in enumerate(signs):
                matrix += (weights[i]*abs(values[i])/ua*(1-acb(0,sign*tilts[i])))*matrices[i]
            tests.append((name,matrix,arb(g),arb(h)))
            bonuses[name] = arb(bonus)
            scores.append((name,matrix,arb(low),L))
    factor = arb("1.75")
    bounds = control["integrated_bounds"](data,tests,bonuses,{},2048,arb(0),budget_factor=factor)
    norm = bounds["norm"]
    assert norm < arb(".016")
    crude = [norm*f for f in lens]
    assert all(v < 1 for v in crude)
    imag = {i:bounds[key+"1"].max(bounds[key+"-1"]) + norm**2*lens[i]/2
            for i,key in ((0,"us"),(1,"cb"),(3,"cs"))}
    assert imag[0] < arb(".04") and imag[1] < arb(".006")
    ub0,ub1 = arb(".003677"),arb(".003787")
    c0,c1 = (1-ub1**2).sqrt(),(1-ub0**2).sqrt()
    mod = ((arb(".224955")*c0,arb(".225065")*c1),
           (arb(".041775")*c0,arb(".041885")*c1),(ub0,ub1))
    q0 = (values[0]*values[1]*values[2].conjugate()*values[3].conjugate()).arg()
    lo,hi = arb("1.141"),arb("1.153")
    offset = sum(w*(abs(values[i])-mod[i][0])/ua for i,w in enumerate(weights[:2])) + weights[2]*(q0-lo)/up
    assert 0 < q0 < arb.pi()/2 and q0+sum(v.asin() for v in crude) < arb.pi()
    assert q0-sum(v.asin() for v in crude) > -arb.pi()
    worst = arb(-100)
    for name,matrix,low,high in scores:
        ratio = [(a/abs(values[i]),b/abs(values[i])) for i,(a,b) in enumerate(mod)]
        ratio.append((1-high*lens[3],1+high*lens[3]))
        angles = {i:(imag[i].min(high*lens[i])/ratio[i][0]).asin() for i in (0,1,3)}
        correction = (mod[0][1]*mod[1][1]*ub1/((1-ub1**2)*abs(values[3])*ratio[3][0])).asin()
        angles[2] = sum(angles.values()) + abs(q0-lo).max(abs(q0-hi)) + correction
        error = sum(abs(1-a).max(abs(1-b))*angles[i]+b*angles[i]**3/6 for i,(a,b) in enumerate(ratio))
        upper = (offset+bounds[name]-bonuses[name]*low+control["nuclear"](matrix)*high**2/2
                 +weights[2]*error/up)
        worst = worst.max(upper.upper())
    print("Center-free budget factor excluded",factor.str(5),"maximum score upper",worst.upper().str(12))
    assert worst<0


def example(data):
    labels = (11415,11430,11432,11389,11448)
    counts = (20547934198306,118869088944549,850933240432901,4716893508610,4932842915634)
    assert sum(counts)==10**15 and 11431 not in labels
    point,scale,eps = data[:3]
    pts = [point(n) for n in labels]
    times = [arb(t)/10**15 for t in counts]
    g = sum(t*x for t,(_,x,_) in zip(times,pts))
    z = g+scale*sum(t*c for t,(_,_,c) in zip(times,pts))
    factor = abs(g/eps).max(abs(z/eps))
    assert factor < arb("2.04")
    blocks = [free["study"]["native"]["plane_response"](h) for h,_,_ in pts]
    half = arb(2)**-29
    largest = arb(0)
    for signs in product((-1,1),repeat=4):
        cube = [v+s*half+arb(0,half) for v,s in zip(free["box"]["centers"],signs)]
        steps = {n:free["ordered"]["evolution"](block,t,cube) for n,block,t in zip(labels,blocks,times)}
        obs = free["ordered"]["readings"](steps,labels[::-1])
        for v,m,w in zip(obs,free["flavor"]["shell"]["centers"],[arb(".000055")]*3+[arb(".006")]):
            largest = largest.max(abs(v-m)/w)
    assert largest<1
    print("Explicit clock budget factor",factor.str(15),"CKM maximum normalized error",largest.upper().str(12))


if __name__ == "__main__":
    data = control["background"]()
    example(data)
    exclude(data)

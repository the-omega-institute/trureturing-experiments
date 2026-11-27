"""逐边可选二阶酉门时，三壳弱端点同时归零的最小平方弦长上限。"""
from itertools import combinations, product
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx


EDGES = ((6,11),(7,12),(8,18),(0,2),(1,3),(4,5))
ANGLES = ("-4.04847308836726","3.87109375","-2.39499703433442","1.37533846482333",
          "-.77821029083598","-1.82024283976827","3.0732180986634","-.28125",
          "2.54124617968621","2.6796875","-1.12109375",".49741970920459",
          "-.3962126605614","-.90234375","-.29296875",".4453125",
          "-.3041420781695","-1.29549312700394","2.06640625","-.4921875",
          "1.3203125","2.4296875","1.83203125","-1.91608546355491")
ACTIVE = (23,3,4,11,17,2,0,8,16,12,5,6)


def adj(a):
    return a.conjugate().transpose()


def experiment(root,fifth,columns,levels):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    rank = lib("three-shell-rank-one-obstruction")
    _,cap = rank["experiment"](root,fifth,columns,levels)
    native = lib("native-rh-ray-pvm")
    geo = lib("three-ray-chiral-boundary")["geo"]
    dist = {(i,j): sum((u-v)**2 for u,v in zip(geo["sites"][i],geo["sites"][j],strict=True))
            for i,j in combinations(range(19),2)}
    assert {e for e in dist if dist[e] < 32} == set(EDGES[:3])
    assert {e for e in dist if dist[e] <=32} == set(EDGES)
    source = fifth*columns
    bases,weights = [],[]
    for n in range(2_306_888_100_092,2_306_888_100_095):
        _,_,_,basis,_,w,_,_ = lib("fifth-seventh-shortest-site-gates")["geometry"](root,fifth,columns,levels,number=n)
        bases.append(adj(basis))
        weights.append(w)
    ctx.prec = 256
    fixed = sorted(set(range(19))-{i for e in EDGES[:3] for i in e})
    scores = []
    cert = root/"ch39-short-pairs-three-weak-separators.tsv"
    selections = product(*rank["CANDIDATES"])
    for line,indices in zip(cert.read_text().splitlines(),selections,strict=True):
        raw = list(map(int,line.split()))
        assert tuple(raw[:3]) == indices and len(raw) == 18
        numbers = raw[6:]
        y = acb_mat([[acb(numbers[4*i+2*j],numbers[4*i+2*j+1])/10000 for j in range(2)] for i in range(3)])
        rows = []
        for b,i,pivot in zip(bases,indices,raw[3:6],strict=True):
            phase = b[i,pivot].conjugate()/abs(b[i,pivot])
            rows.append([phase*b[i,j] for j in range(19)])
        s = acb_mat(rows)
        def pieces(indices):
            return (acb_mat([[s[h,i] for i in indices] for h in range(3)]),
                    acb_mat([[source[i,j] for j in range(2)] for i in indices]))
        rest,fixed_columns = pieces(fixed)
        a = adj(y)*rest*fixed_columns
        margin = (a[0,0]+a[1,1]).real
        for edge in EDGES[:3]:
            se,ce = pieces(edge)
            g = adj(se)*y*adj(ce)
            margin -= (sum(abs(z)**2 for z in g.entries())+2*abs(g.det())).sqrt()
        ynorm = sum(abs(z)**2 for z in y.entries()).sqrt()
        scores.append((margin/ynorm).lower())
    print("三条最短边在十八种低线选择下的投影范数一致下界",min(scores).str(13),flush=True)
    assert min(scores) > arb(".014") > cap

    gate = lib("two-shell-ordered-bond-gates")["transport"]
    layers = tuple((k,(e,)) for k,e in enumerate(EDGES))
    chosen = (0,9,15)
    rows = acb_mat([[base[k,i] for i in range(19)] for base,k in zip(bases,chosen,strict=True)])
    weak = native["weak_series"](56)
    low,heavy,top = levels
    mp,mz = arb("1.22089e19"),arb("91.1876")
    p0,p1 = (mz/mp)**2,(arb("1.176397e14")/mp)**2
    ref = arb("-.074256577194663365545")
    large = lib("full-neutral-free-cutoff")["weak_column"]
    qcenters = []
    for w,k in zip(weights,chosen,strict=True):
        left,right = arb(105),arb(2000)
        for _ in range(48):
            mid = (left+right)/2
            if (ref+w[k]*(weak((mid/mz)**2)-ref)).mid() > 0:
                left = mid
            else:
                right = mid
        qcenters.append((left+right).mid()/2)
    width = arb(".001")
    qs = [q+arb(0,width) for q in qcenters]

    def matrices(point):
        return gate(source,point,parameters=ANGLES,active=ACTIVE,layers=layers)

    def system(point,q):
        mats = [rows*m for m in matrices(point)]
        f,j = [],[]
        for h,energy in enumerate(q):
            for col in range(2):
                goal = ((energy/mp-low)/(heavy-low)).sqrt() if col==0 else arb(0)
                z = mats[0][h,col]-goal
                f.extend(([z.real],[z.imag]))
                j.extend(([m[h,col].real for m in mats[1:]], [m[h,col].imag for m in mats[1:]]))
        return arb_mat(f),arb_mat(j)

    center = [arb(ANGLES[k]) for k in ACTIVE]
    f,j = system(center,qcenters)
    step = j.mid().solve(f.mid()).entries()
    center = [(q-v).mid() for q,v in zip(center,step,strict=True)]
    f,j = system(center,qs)
    pre = j.mid().inv().mid()
    assert not pre.det().contains(0)
    radius = arb(".00005")
    box = [x+arb(0,radius) for x in center]
    _,jbox = system(box,qs)
    eye = arb_mat([[int(i==j) for j in range(12)] for i in range(12)])
    ratio = max((v.abs_upper()/radius+sum(z.abs_upper() for z in r)).upper()
                for v,r in zip((pre*f).entries(),(eye-pre*jbox).tolist(),strict=True))
    print("六条不交边的一轮十二式闭盒自映比",ratio.str(12),flush=True)
    assert ratio < arb(".064")
    operated = matrices(box)[0]
    for k,w,b,q in zip(chosen,weights,bases,qcenters,strict=True):
        p = b*operated
        remaining = [low+sum((levels[j+1]-low)*abs(p[i,j])**2 for j in range(2)) for i in range(19) if i!=k]
        minimum = min(a.lower() for a in remaining)
        assert minimum > 0 and minimum**2 > p0 and sum(w).contains(1)
        tail = sum(w[i]*(large(minimum**2,p0,p1)-ref) for i in range(19) if i!=k)
        faces = [ref+w[k]*(weak(((q+sgn*width)/mz)**2)-ref)+tail for sgn in (-1,1)]
        print("所选谱线、能标中心、相对两面的弱差",k+1,q.str(15),faces[0].str(10),faces[1].str(10),flush=True)
        assert faces[0] > 0 > faces[1]


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    source = run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    experiment(root,*source["verify"](root,return_operator="spectral"))

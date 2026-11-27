"""两套十二短边九式见证在三个原校准壳的中性规范残差与匹配点色动能。"""
from pathlib import Path
from runpy import run_path

from flint import arb


def responses(bands, weights, data, species, gauge, rule):
    p0,p1 = data["p0"],data["p1"]
    log = (p1/p0).log()
    al,aw,ac = data["initial"]
    values = [5*al/3+aw-8*ac/3,aw-ac,ac]
    for (left,right),weight in zip(bands,weights,strict=True):
        assert 0 < left <= right
        for kind,mass,(q,m,c) in species:
            coeff = (q-arb(8)*c/3,m-c,c)
            pair = [gauge["angular_exposure"](kind,mass,end,p0,p1,log,rule) for end in (left,right)]
            exposure = pair[0].union(pair[1])*weight/(4*arb.pi())
            for j,b in enumerate(coeff):
                values[j] -= b*exposure
    return values


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    models = [lib("twelve-independent-gates-nine"),lib("second-twelve-gates-nine")]
    transport = lib("two-shell-ordered-bond-gates")["transport"]
    layers = tuple((j,(edge,)) for j,edge in enumerate(models[0]["PAIRS"]))
    basis_cells = []
    source = fifth*columns
    for model in models:
        box,energies = model["experiment"](root,fifth,columns,levels)
        def mats(point):
            return transport(source,point,parameters=model["PARAMETERS"],active=model["ACTIVE"],layers=layers)
        center = [v.mid() for v in box]
        radii = [(v-c).abs_upper() for v,c in zip(box,center,strict=True)]
        basis_cells.append((mats(center)[0],mats(box)[1:],radii,energies))
    mp = arb("1.22089e19")
    rho = (arb("246.22")/mp)**2/2
    low,heavy,top = levels
    assert top-low < arb(".756")
    gauge = lib("all-field-signed-crossing")
    rule = gauge["angular_rule"](13,39)
    output = []
    for shell,number in enumerate(range(2_306_888_100_092,2_306_888_100_095)):
        data = gauge["joint_system"](root,"1.176397e14",number=number,order=4,last=8)
        signs,masses = lib("rotation-portal-top-energy")["rotation_signs"](root,number)
        species,_ = lib("full-neutral-free-cutoff")["response_rows"](data,masses,signs,rho,arb(".42995401872"))
        _,_,_,basis,_,weights,_,_ = lib("fifth-seventh-shortest-site-gates")["geometry"](root,fifth,columns,levels,number=number)
        left_basis = basis.conjugate().transpose()
        for model,(center,derivatives,radii,energies) in zip(models,basis_cells,strict=True):
            proj = left_basis*center
            dirs = [left_basis*m for m in derivatives]
            bands = []
            for i in range(19):
                if i==model["LINES"][shell]:
                    energy = energies[shell]/mp
                    bands.append((energy.lower()**2,energy.upper()**2))
                    continue
                endpoints = []
                for c in range(2):
                    motion = sum(r*abs(d[i,c]).upper() for r,d in zip(radii,dirs,strict=True))
                    z = abs(proj[i,c])
                    endpoints.append((max(arb(0),(z.lower()-motion).lower()),(z.upper()+motion).upper()))
                a,b = [low+sum(gain*ends[k]**2 for gain,ends in zip((heavy-low,top-low),endpoints,strict=True)) for k in (0,1)]
                bands.append((a.lower()**2,min(top.upper(),b.upper())**2))
            vals = responses(bands,weights,data,species,gauge,rule)
            # 算子范数移动 δ 时，每个固定单位本征行上的对角至多移动 δ。
            delta = arb("1e-14")
            wide = [(max(low.lower(),(x.sqrt()-delta).lower())**2,
                     min(top.upper(),(y.sqrt()+delta).upper())**2) for x,y in bands]
            nearby = responses(wide,weights,data,species,gauge,rule)
            print("算子范数闭邻域",delta.str(5),[(v.lower().str(9),v.upper().str(9)) for v in nearby],flush=True)
            assert nearby[0]>8844 and nearby[1]>645 and nearby[2]<-2411
            assert vals[0]>0 and vals[1]>0 and vals[2]<0
            print(number,model["LINES"],"光子减 8/3 色、混合减色、色逆耦合",
                  [(v.lower().str(13),v.upper().str(13)) for v in vals],flush=True)
            output.append(vals)
    for n in range(3):
        first,second = output[2*n:2*n+2]
        assert all(not a.overlaps(b) for a,b in zip(first,second,strict=True))
    return output


def verify(root):
    joint = run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    return experiment(root,*joint["verify"](root,return_operator="spectral"))


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)

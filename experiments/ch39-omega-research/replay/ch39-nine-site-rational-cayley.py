"""从有限整数表构造严格有理的九站酉，直接读取十三项弱差。"""
import json
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx, fmpq, fmpq_mat


def rational_gate(root):
    data = json.loads((root/"ch39-nine-site-rational-cayley.json").read_text())
    sites = data["sites"]
    assert len(sites)==len(set(sites))==9 and min(sites)>=0 and max(sites)<19
    assert data["bits"]==54 and data["quarter_turns"] in range(4)
    denominator = 2**data["bits"]
    h = fmpq_mat(18,18)
    def set_entry(i,j,a,b):
        h[i,j] = h[i+9,j+9] = fmpq(a,denominator)
        h[i+9,j] = fmpq(b,denominator)
        h[i,j+9] = -h[i+9,j]
    assert len(data["upper"]) == 9
    for i,row in enumerate(data["upper"]):
        assert len(row)==9-i
        for j,(a,b) in enumerate(row,start=i):
            assert type(a) is int and type(b) is int
            assert i!=j or a==0
            set_entry(i,j,a,b)
            set_entry(j,i,-a,b)
    assert h.transpose()==-h
    identity = fmpq_mat([[int(i==j) for j in range(18)] for i in range(18)])
    u = (identity-h)*(identity+h).inv()
    assert u.transpose()*u==identity
    phase = acb(0,1)**data["quarter_turns"]
    gate = acb_mat([[phase*acb(arb(u[i,j]),arb(u[i+9,j])) for j in range(9)] for i in range(9)])
    return sites, gate


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    ctx.prec = 320
    sites, local_gate = rational_gate(root)
    assert set(sites)=={1,3,4,5,8,9,10,11,14}
    gate = acb_mat([[int(i==j) for j in range(19)] for i in range(19)])
    for i, site in enumerate(sites):
        for j, target in enumerate(sites):
            gate[site,target] = local_gate[i,j]

    graph = lib("three-ray-chiral-boundary")
    native = lib("native-rh-ray-pvm")
    gram = native["gram_cholesky"](root)[3]
    shell_rows, weights = [], []
    for shell in range(7,20):
        _, basis = lib("adjacent-pvm-transport")["spectrum"](graph,6,2_306_888_100_085+shell)
        shell_rows.append(basis.conjugate().transpose()*gate)
        weights.append([z.real for z in native["readout"](basis,range(19),[14,8,7,6],gram,3)])

    probes = acb_mat([matrix.tolist()[0] for matrix in shell_rows])
    fifth, columns, levels, first_lines = lib("fifth-full-positive-joint")["verify"](
        root, return_operator="spectral", spectral_probes=probes)
    floor, middle, top = levels
    mp, mz = arb("1.22089e19"), arb("91.1876")
    p0, p1 = (mz/mp)**2, (arb("1.176397e14")/mp)**2
    low = native["weak_series"](56)
    high = lib("full-neutral-free-cutoff")["weak_column"]
    ref = arb("-.074256577194663365546").union(arb("-.074256577194663365544"))
    mass_weights = (middle-floor,top-floor)
    weighted = columns*acb_mat([[mass_weights[0].sqrt(),0],[0,mass_weights[1].sqrt()]])
    bounds = []
    for shell,(matrix,w,first) in enumerate(zip(shell_rows,weights,first_lines.tolist(),strict=True),7):
        upper_lines = ((matrix*fifth)*weighted).tolist()[1:]
        assert all(z.is_finite() for row in upper_lines for z in row)
        minimum = min((floor+sum(max(arb(0),abs(z).lower())**2 for z in row)).lower()
                      for row in upper_lines)
        assert minimum>0 and minimum**2>p0
        tail = (high(minimum**2,p0,p1)-ref).abs_upper()*sum(abs(z).upper() for z in w[1:])
        q = mp*(floor+sum(k*abs(z)**2 for k,z in zip(mass_weights,first,strict=True)))
        d = ref+w[0]*(low((q/mz)**2)-ref)+arb(0,tail)
        assert d.is_finite()
        print("实际壳、第一线截止 GeV 与有理矩阵弱差",shell,q.str(10),d.str(12),flush=True)
        bounds.append(abs(d).upper())
    print("十三壳最大弱差绝对上界",max(bounds).str(13),flush=True)
    assert max(bounds)<arb(".00076")
    return max(bounds)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)

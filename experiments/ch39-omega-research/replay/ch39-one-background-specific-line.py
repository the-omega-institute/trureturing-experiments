"""两组实际背景只在一条线上独立分配权重，精确匹配双引力及正交汇。"""
import json
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

ROOT = Path(__file__).resolve().parent
LABELS = (
    "000bb10000b1b001bab","2227600b811b300b101","b2311b0bb12b410b000",
    "14077006400b100b000","11176008601b200b000","1127700b912b300b000",
    "0bbbb0b88690845b545","2bb890752bbbbb5b000","b78344310ab10820bbb",
    "0bb98062100b100b000","12311b11048bb41b000","18a45b442bbbbb4b000",
    "14011b020b0b000b000","b12019bbbb509bb0bbb","10000b032b0b1b0b000",
    "10000bbbbb0b0bb0000","10000bb0bb0b0b0b000","0bbbb03000b0b0b0bbb",
    "10000bb0bb0b0b0b000",
)
# (背景或 None=共用，线编号，额外格点)。每个自由份额来自该线的登记基点。
EXTRA = ((None,17,11),(None,285,2),(None,285,11),(None,302,11),(None,323,10),
         (0,0,11),(1,0,6),(1,0,7))


def load(n):
    return run_path(str(ROOT/f"ch39-{n}.py"))


def solve_layout(labels, extra, independent, minima, shell_codes, heavy_codes, fit_energy,
                 *, reference_gev=10**16, scale_millionths=None, weight_codes=None, positive_codes=None, edge_tables=None):
    source = load("shared-twelve-node-gauge-gap")
    single,lib,backgrounds = load("two-zero-cutoff-gravity")["gravity_fields"](
        ROOT,heavy_codes,tight_neutral=True)
    model,_,_,initial,_ = lib["chiral"]["source"]["mean"]["inputs"]()
    gauge = load("signed-gravity-unification")
    shell = load("native-rh-ray-pvm")["shell_weights"]
    fields = gauge["gauge_weights"](model)
    beta = [b/(2*arb.pi()) for b in model["beta"]]
    logscale = (arb(reference_gev)/model["mz"]).log()
    bare = [v-b*logscale for v,b in zip(initial,beta,strict=True)]
    count = len(shell_codes)
    scale_millionths = [7000]*count if scale_millionths is None else scale_millionths
    weight_codes = shell_codes if weight_codes is None else weight_codes
    positive_codes = shell_codes if positive_codes is None else positive_codes
    assert len(weight_codes)==len(positive_codes)==count
    independent = {independent} if isinstance(independent,int) else set(independent)
    assert len(scale_millionths)==count and all(v>0 for v in scale_millionths)
    assert len(labels)==count and len(extra)+int(fit_energy)==4*count
    assert all(len(row)==361 for row in labels)
    for i,row in enumerate(zip(*labels,strict=True)):
        if i not in independent:
            assert all(a is row[0] or a==row[0] for a in row[1:])
    assert all(who is None or (i in independent and who in range(count)) for who,i,k in extra)
    bases,changes = [],[]
    for n,(code,bg) in enumerate(zip(shell_codes,backgrounds,strict=True)):
        if edge_tables is None:
            edges = source["paired_lines"](ROOT,shell,code,weight_code=weight_codes[n],positive_code=positive_codes[n])
        else:
            assert len(edge_tables)==count
            edges = edge_tables[n]
            assert len(edges)==361 and all(0<a and (a<b or a==b) and b<=1 for a,b,w in edges)
        scale = arb(scale_millionths[n])/10**6
        def response(i,k):
            a,b,w = edges[i]
            assert 0 <= k <= 1
            x = scale*((1-k)*a+k*b)
            _,K,W,g1,g2,g3 = gauge["threshold_column"](bg,single["moments"],lib["expint"],fields,x)
            return arb_mat([[w*z] for z in (K,W,g1-g3,g2-g3,g3)])
        chosen = [response(i,k) for i,k in enumerate(labels[n])]
        base = arb_mat([[-1],[-arb("1e-122")],[bare[0]-bare[2]],[bare[1]-bare[2]],[bare[2]]])
        for row in chosen:
            base += row
        bases.append(base.entries())
        changes.append([(response(i,k)-chosen[i]).entries() if who in (None,n) else [arb(0)]*5
                        for who,i,k in extra])
    columns = [[v for group in changes for v in group[j][:4]] for j in range(len(extra))]
    if fit_energy:
        columns.append([0,0,beta[2]-beta[0],beta[2]-beta[1]]*count)
    matrix = arb_mat(columns).transpose()
    rhs = arb_mat([[-v] for base in bases for v in base[:4]])
    def midpoint(m):
        return arb_mat([[m[i,j].mid() for j in range(m.ncols())] for i in range(m.nrows())])
    pre = midpoint(midpoint(matrix).inv())
    seed = midpoint(pre*midpoint(rhs))
    normalized = pre*matrix
    assert not normalized.det().contains(0)
    solution = (seed+normalized.solve(pre*(rhs-matrix*seed))).entries()
    shares = solution[:len(extra)]
    log_change = solution[-1] if fit_energy else arb(0)
    energy = arb(reference_gev)*log_change.exp()
    assert model["mz"] < energy < arb("1.22089e19")
    assert all(0 < s < 1 for s in shares)
    for n in range(count):
        unused = {i:arb(1) for _,i,_ in extra}
        for s,(who,i,_) in zip(shares,extra,strict=True):
            if who in (None,n):
                unused[i] -= s
        assert all(v > 0 for v in unused.values())
        common = bases[n][-1]+sum(s*row[-1] for s,row in zip(shares,changes[n],strict=True))
        minimum = common+beta[2]*logscale
        common -= beta[2]*log_change
        assert beta[2]<beta[1]<0<beta[0] and common>0 and minimum>arb(minima[n])
        print("背景、共同逆耦合、从电弱到交汇能标的最小逆耦合",n,common.str(13),minimum.str(13))
    print("内部格点份额",[x.str(11) for x in shares])
    return shares,energy


def certify(labels, extra, independent, minima):
    source = load("shared-twelve-node-gauge-gap")
    shares,_ = solve_layout(labels,extra,independent,minima,source["SHELL"],source["HEAVY"],False)
    return shares


def layouts():
    labels = [arb(int(d,12))/11 for d in "".join(LABELS)]
    yield [labels,labels],[(n,i,arb(k)/11) for n,i,k in EXTRA],0,("1.535","1.231")
    layout = json.loads((ROOT/"ch39-background-specific-line-323.json").read_text())
    denominator = layout["denominator"]
    labels = [arb(v)/denominator for row in layout["common_base_numerators"] for v in row]
    second = list(labels)
    for i,v in layout["second_base_changes"]:
        second[i] = arb(v)/denominator
    extra = [(n,i,arb(v)/denominator) for n,i,v in layout["extra_numerators"]]
    yield [labels,second],extra,323,("1.528","1.056")


def verify():
    for layout in layouts():
        certify(*layout)


if __name__ == "__main__":
    verify()

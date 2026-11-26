"""三条腿各取原场或一次邻边输运时的汤川质量与混合。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


edge = run_path(str(Path(__file__).with_name("ch39-single-edge-yukawa.py")))
shell = edge["shell"]
adjoint = edge["filters"]["local"]["adjoint"]
masks = tuple(product(range(2), repeat=3))


def sector(left, right, higgs, conjugate):
    rows = [[] for _ in range(3)]
    for k in range(3):
        for i in range(3):
            j = (-i + (k if conjugate else -k)) % 3
            rows[k].append([
                sum(left[a][x,i] * right[b][x,j] * higgs[c][x,k] for x in range(19))
                for a,b,c in masks
            ])
    return [acb_mat(entries) for entries in rows]


def tensors(number):
    height = acb.zeta_zero(number).imag
    neutral = shell["laplacian"](0, height)
    center = sum(neutral[i,i].real for i in range(19)) / 18
    fields = {}
    transports = {}
    heavy_minima = []
    for q in (1, -4, 2, -3, 6, 3):
        lap = shell["laplacian"](q, height)
        transport = acb_mat([[-lap[i,j]/8 if i != j else 0 for j in range(19)] for i in range(19)])
        modes = edge["mode_columns"](edge["filters"]["spectra"](q, height, center))
        # 固定低模相位，使同一见证在解压后的重放中也有相同坐标。
        for j in range(3):
            anchor = next(modes[i,j] for i in range(19) if abs(modes[i,j]) > arb(".01"))
            phase = anchor.conjugate()/abs(anchor)
            for i in range(19):
                modes[i,j] *= phase
        fields[q] = (modes, transport * modes)
        transports[q] = transport
        if q in (1,-4):
            roots = lap.eig(algorithm="rump")
            levels = [abs(r.real*(r.real-center)) for r in roots]
            assert all(v > 9 or v < 9 for v in levels)
            heavy = [v for v in levels if v > 9]
            assert len(heavy) == 16
            unit = arb("2.435e18")/724**2
            heavy_minima.append(min(v.lower() for v in heavy)*unit*(32 if q==-4 else 1))
    hbar = tuple(block.conjugate() for block in fields[3])
    blocks = list(zip(sector(fields[1], fields[-4], fields[3], False),
                    sector(fields[1], fields[2], hbar, True),
                    sector(fields[-3], fields[6], hbar, True)))
    return blocks, full_up_gram([transports[q] for q in (1,-4,3)]), higgs_up_grams(transports, fields), heavy_minima


def full_up_gram(transports):
    identity = acb_mat([[int(i==j) for j in range(19)] for i in range(19)])
    kernels = []
    for transport in transports:
        ops = (identity, transport)
        kernels.append({(i,j): a.conjugate()*b.transpose() for i,a in enumerate(ops) for j,b in enumerate(ops)})
    gram = acb_mat(8,8)
    for i,a in enumerate(masks):
        for j,b in enumerate(masks):
            left,right,higgs = [ker[x,y] for ker,x,y in zip(kernels,a,b)]
            gram[i,j] = sum(left[s,t]*right[s,t]*higgs[s,t] for s in range(19) for t in range(19))
    return gram


def higgs_up_grams(transports, fields):
    higgs = fields[3]
    vq, vu = fields[1][0], fields[-4][0]
    identity = acb_mat([[int(i==j) for j in range(19)] for i in range(19)])
    left, right = [(identity, transports[q]) for q in (1,-4)]
    grams = []
    heavy_left = identity - vq.conjugate()*vq.transpose()
    heavy_right = identity - vu*adjoint(vu)
    for k in range(3):
        columns, heavy_columns = [], []
        for a,b,c in masks:
            h = higgs[c]
            column = left[a].transpose()*acb_mat([[h[x,k]*right[b][x,j] for j in range(19)] for x in range(19)])
            heavy = heavy_left*column*heavy_right
            columns.append(list(column.entries()))
            heavy_columns.append(list(heavy.entries()))
        coeff, heavy_coeff = [acb_mat(rows).transpose() for rows in (columns, heavy_columns)]
        grams.append((adjoint(coeff)*coeff, adjoint(heavy_coeff)*heavy_coeff))
    return grams


def gram_pair(block):
    u,d,e = [adjoint(c)*c for c in block]
    return u, (arb("61.36519")**2*d + arb("84.98147")**2*e) / 2


WITNESS = (
    ("8.545351061704176", "-48.47372271439794"),
    ("2.2078900661717036", "165.2152361507395"),
    ("9.60539599947788", "-82.330234253038"),
    ("-822.4707615614616", "3317.586949639615"),
    ("0.16611420736276372", "-139.0101528268871"),
    ("-578.8448874434607", "2943.6299896262976"),
    ("848.4310720697915", "-3172.7679246028706"),
    ("-896.7912251591839", "4432.855036558063"),
)
HIGGS = (
    ("0.9999999915763772", "0.0"),
    ("8.214850466378228e-05", "9.957579069062194e-05"),
    ("-1.1792222338125794e-05", "-6.668894941105885e-06"),
)


def column(values):
    return acb_mat([[acb(arb(x),arb(y))] for x,y in values])


def positive(matrix):
    for size in range(1, matrix.nrows()+1):
        minor = acb_mat([[matrix[i,j] for j in range(size)] for i in range(size)])
        assert minor.det().real > 0


def mass_matrices(blocks, w, z):
    matrices = [acb_mat(3,3) for _ in range(3)]
    for k, block in enumerate(blocks):
        for f, coeff in enumerate(block):
            entries = coeff * w
            for i in range(3):
                j = (-i + (-k if f == 0 else k)) % 3
                amplitude = z[k,0] if f == 0 else z[k,0].conjugate()
                matrices[f][i,j] = amplitude*entries[i,0]
    return matrices


def verify():
    refs = run_path(str(Path(__file__).with_name("ch39-shell-flavor-matching-windows.py")))
    mixing = run_path(str(Path(__file__).with_name("ch39-finite-orthogonal-ckm.py")))
    ctx.prec = 448
    w, z = column(WITNESS), column(HIGGS)
    z = z/(adjoint(z)*z)[0,0].real.sqrt()
    all_deficits = []
    for n, cap, full_cap, heavy_cap in ((22027,".000000913",".000127",".000167"),
                                        (28325,".000000890",".000126",".000166")):
        blocks, full, physical, minima = tensors(n)
        pairs = [gram_pair(b) for b in blocks]
        deficits = [d-u for u,d in pairs]
        for k,l in product(range(3), repeat=2):
            if k != l:
                positive(deficits[k] + arb(".000061")*deficits[l])
        all_deficits.append(deficits)
        for (u,d), (both,heavy) in zip(pairs,physical):
            positive(arb(cap)*full - u + arb(".015")*(d-u))
            positive(arb(full_cap)*both - u + arb(".088")*(d-u))
            positive(arb(heavy_cap)*heavy - u + arb(".085")*(d-u))
        assert refs["up"][-1]*arb(".999") > arb(".86")
        print(n, "给定真空方向的完整、重重上型矩阵谱范数下界",
              (arb(".86")/(19*arb(full_cap)).sqrt()).str(9),
              (arb(".86")/(16*arb(heavy_cap)).sqrt()).str(9))
        q2, u2 = [v**2 for v in minima]
        heavy_norm2 = arb(".86")**2/arb(heavy_cap)
        quadratic = (q2+u2-q2*u2/(q2+u2))*heavy_norm2
        quartic = heavy_norm2**2/64
        assert quadratic > arb("3.45e34") and quartic > 306000
        print(n, "最轻重双重态与单态 GeV、四次质量迹的二次及四次系数下界",
              *(v.str(9) for v in (*minima,quadratic,quartic)))
        floor = arb(".86")/(19*arb(cap)).sqrt()
        print(n,"完整上型顶角展开的谱范数下界",floor.str(9))
        matrices = mass_matrices(blocks,w,z)
        target = [v for name in ("up","down","leptons") for v in refs[name]]
        masses = [s for y in matrices for s in sorted(
            (r.real.sqrt() for r in (y*adjoint(y)).eig(algorithm="rump")), key=float)]
        errors = [abs(v/t-1) for v,t in zip(masses,target)]
        angles = mixing["observables"](*matrices[:2])[4:]
        if n == 22027:
            expected = [arb(x).sin() for x in (".22704",".04275",".00377")] + [arb("1.139")]
            assert all(e < arb("3.3e-8") for e in errors)
            assert all(abs(v-ref) < arb("1e-8") for v,ref in zip(angles[:4], expected))
            assert arb(".00003206") < angles[4] < arb(".00003207")
        print(n,"九个质量相对误差中的最大值",max(err.upper() for err in errors).str(9))
        print(n,"混合角正弦、CP 相位与 J",*(v.str(10) for v in angles))
        by_char = []
        for block in blocks:
            norms = [(adjoint(b*w)*(b*w))[0,0].real.sqrt() for b in block]
            by_char.append((norms[0]/norms[1], norms[0]/norms[2]))
        if n == 28325:
            assert all(e < arb("47.366") for _,e in by_char)
        print(n,"固定八系数任意 Higgs 的两范数比上限",
              *(max(v.upper() for v in row).str(10) for row in zip(*by_char)))
    for k,l in product(range(3), repeat=2):
        if k != l:
            positive(all_deficits[0][k] + all_deficits[1][l])


if __name__ == "__main__":
    verify()

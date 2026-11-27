"""固定第五壳谱，一个三循环等变酉输运可以在第七至十二壳同时实现六个零弱差。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


Q = ("522.398135900", "534.208610459", "544.705683199",
     "544.755583188", "533.193826337", "516.690391112")


def adj(m):
    return m.conjugate().transpose()


def squared(m):
    return sum(abs(z)**2 for z in m.entries())


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    graph = lib("three-ray-chiral-boundary")
    native = lib("native-rh-ray-pvm")
    characters, rot = lib("cyclic-two-characters-nine")["characters"](graph["geo"]["sites"])
    trivial = characters[0]
    # 坐标循环保持球面内积与行列式，带钟边不反向，软辐边组成一个轨道。
    edges = {(a,b):(d,z) for a,b,d,z in graph["geo"]["edges"]}
    assert rot[graph["soft"]["pole"]] == graph["soft"]["pole"]
    orbit = set(graph["soft"]["star"]["orbit"])
    assert {rot[k] for k in orbit} == orbit
    for (a,b),(d,z) in edges.items():
        i,j = rot[a],rot[b]
        e,w = edges[tuple(sorted((i,j)))]
        assert e == d and (d not in (724,804) or i<j)
        assert (w if i<j else w.conjugate()).overlaps(z)

    gram = native["gram_cholesky"](root)
    bases = []
    for j in range(6):
        _, v = lib("adjacent-pvm-transport")["spectrum"](graph, 6, 2_306_888_100_092+j)
        weights = [w.real for w in native["readout"](v, range(19), [14,8,7,6], gram[3], 3)]
        bases.append((adj(v),weights))
    ctx.prec = 256
    rows = acb_mat([[basis[0,k] for k in range(19)] for basis,_ in bases])
    b = rows*trivial
    dual = adj(b)*(b*adj(b)).inv()
    null = acb_mat([[int(i==j) for j in range(7)] for i in range(7)])-dual*b
    index = max(range(7), key=lambda k: null[k,k].real.lower())
    n = acb_mat([[null[k,index]] for k in range(7)])
    n = n/squared(n).sqrt()
    floor,heavy,top = levels
    weighted = (fifth*columns)*acb_mat([[(heavy-floor).sqrt(),0],[0,(top-floor).sqrt()]])
    t = adj(trivial)*weighted
    gram2 = adj(t)*t
    trace = gram2[0,0].real+gram2[1,1].real
    det = sum(abs(t[i,0]*t[j,1]-t[i,1]*t[j,0])**2 for i in range(7) for j in range(i))
    high = (trace+(trace**2-4*det).sqrt())/2
    low = det/high

    planck, mz = arb("1.22089e19"), arb("91.1876")
    width = arb(".00005")
    centers = list(map(arb,Q))
    y = acb_mat([[((v+arb(0,width.upper()))/planck-floor).sqrt()] for v in centers])
    lift = dual*y
    m = squared(lift)
    assert low < m < high
    u = lift/m.sqrt()
    coupling = ((high-m)*(m-low)).sqrt()
    ref = arb("-.074256577194663365546").union(arb("-.074256577194663365544"))
    weak = native["weak_series"](56)
    large = lib("full-neutral-free-cutoff")["weak_column"]
    p0, p1 = (mz/planck)**2,(arb("1.176397e14")/planck)**2
    print("两非零本征值与特解范数平方",*[v.str(12) for v in (low,m,high)],flush=True)
    for j,(basis,weights) in enumerate(bases):
        others=[]
        for i in range(19):
            row=acb_mat([[basis[i,k] for k in range(19)]])
            f=row*trivial
            # 原生单谱行只属一个字符，用非零分量辨认它，避免由区间含零断定精确零。
            memberships=[k for k,E in enumerate(characters) if any(abs(z)>arb(".01") for z in (row*E).entries())]
            assert len(memberships)==1
            if i==0:
                assert memberships==[0]
                continue
            if memberships==[0]:
                a, c = (f*u)[0,0], (f*n)[0,0]
                amp = floor + m*abs(a)**2 + (trace-m)*abs(c)**2 + 2*coupling*(a*c.conjugate()).real
            else:
                amp = floor+squared(row*weighted)
            others.append(amp.lower())
        minimum = min(others)
        assert minimum**2 > p0 and minimum > 0,(j,minimum)
        tail = (large(minimum**2,p0,p1)-ref).abs_upper()*sum(abs(w).upper() for w in weights[1:])
        faces=[ref+weights[0]*(weak(((centers[j]+sign*width)/mz)**2)-ref)+arb(0,tail)
               for sign in (-1,1)]
        print("壳、其他线最小截止下界及弱差的两个面",j+7,minimum.str(12),faces[0].lower().str(12),faces[1].upper().str(12),flush=True)
        assert faces[0]>0>faces[1]
    # 基 (u,n) 上取 Hermitian 块 [[m,coupling],[coupling,trace-m]]。它与源的平凡字符
    # 块同谱，且 Bn=0、Bu=y/sqrt(m)，故六条低线正好具有各自指定的截止。
    # 非平凡字符块取恒等。读数连续且六对相反面变号，由 Poincaré–Miranda 得共同零点。


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    src = run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    experiment(root, *src["verify"](root, return_operator="spectral"))

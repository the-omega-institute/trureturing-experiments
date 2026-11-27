"""两个实际零点的二次非交换 Gram 多项式源，以及整盒两相位存在性证书。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat


def selected_rows(matrix, indices):
    return arb_mat([[matrix[i,j] for j in range(matrix.ncols())] for i in indices])


def mixed_words(up, down):
    return (down, up*down+down*up, acb(0,1)*(up*down-down*up), down*down)


def certify_box(orbit, model, up_map, down_map):
    adj = model["adj"]
    values, frame = model["eig"](adj(up_map)*up_map)
    up = acb_mat([[v if i==j else 0 for j in range(3)] for i,v in enumerate(values)])
    down = adj(frame)*(adj(down_map)*down_map)*frame
    pairs = ((0,1),(0,2),(1,2))
    assert all(abs(down[i,j]) > 0 for i,j in pairs)

    def coordinates(matrix):
        z = [matrix[i,j]/down[i,j] for i,j in pairs]
        return arb_mat([[v.real] for v in z]+[[v.imag] for v in z])

    columns = [coordinates(word) for word in mixed_words(up,down)]
    features = arb_mat([[column[i,0] for column in columns] for i in range(6)])
    leading, remaining = (0,2,3,4), (1,5)
    minor = selected_rows(features,leading)
    assert abs(minor.det()) > 0
    projection = selected_rows(features,remaining)*minor.inv()

    def residue(matrix):
        vector = coordinates(matrix)
        return selected_rows(vector,remaining)-projection*selected_rows(vector,leading)

    direct,_ = orbit["orbit"](model,up_map,down_map)
    parts = [adj(frame)*part for part in direct]
    constant = residue(sum((a*adj(a) for a in parts),acb_mat(3,3)))
    terms = []
    for i,j in pairs:
        cross = parts[i]*adj(parts[j])
        cosine = residue(cross+adj(cross))
        sine = residue(acb(0,1)*(cross-adj(cross)))
        terms.append((int(i==0)-int(j==0),int(i==1)-int(j==1),cosine,sine))

    def residual_jacobian(x,y):
        value = constant
        dx,dy = arb_mat(2,1),arb_mat(2,1)
        for a,b,cosine,sine in terms:
            t = a*x+b*y
            value = value + cosine*t.cos() + sine*t.sin()
            slope = -cosine*t.sin()+sine*t.cos()
            dx = dx + a*slope
            dy = dy + b*slope
        return value,arb_mat([[dx[0,0],dy[0,0]],[dx[1,0],dy[1,0]]])

    x,y = arb("1.75934538"),arb("4.34780044")
    half = arb(".05")
    offset = arb(0,half.upper())
    f,center_jacobian = residual_jacobian(x,y)
    midpoint = arb_mat([[center_jacobian[i,j].mid() for j in range(2)] for i in range(2)])
    assert abs(midpoint.det())>0
    inverse = midpoint.inv()
    correction = max(abs(v).upper() for v in (inverse*f).entries())
    _,jacobian = residual_jacobian(x+offset,y+offset)
    derivative = arb_mat([[1,0],[0,1]])-inverse*jacobian
    lipschitz = max(sum(abs(derivative[i,j]).upper() for j in range(2)) for i in range(2))
    assert correction < arb(".0259") and lipschitz < arb(".417")
    assert correction + half*lipschitz < half
    print("整盒根存在，中心位移上界和导数行范数上界",correction.str(12),lipschitz.str(12))


def at_zero(model,q,n):
    adj = model["adj"]
    laplacian = model["shell"]["laplacian"](q,acb.zeta_zero(n).imag)
    rows = []
    for frame in model["frames"]:
        levels, basis = model["eig"](adj(frame)*laplacian*frame)
        assert levels[0] < arb("8.685") and levels[1] > arb("9.6")
        first = acb_mat([[basis[i,0]] for i in range(basis.nrows())])
        row = adj(frame*first)*model["sources"]/(1+(levels[0]/9)**512)
        rows.append(row.entries())
    return acb_mat(rows)


def rational_covariance(u, d):
    identity = acb_mat([[int(i==j) for j in range(3)] for i in range(3)])
    words = (identity,u,u*u,*mixed_words(u,d))
    numerators = (15570864405642,-625228858823767,-12873069640991703,
                  -10533333541133613,16612062805942317,-90739907284335,47162635274051528)
    return sum((word*(arb(n)/10**12) for word,n in zip(words,numerators,strict=True)),
               acb_mat(3,3))


def fixed_polynomial(model):
    adj = model["adj"]
    sites = [[at_zero(model,charge,index+offset) for offset in (-1,0,1)]
             for charge,index in ((-4,2345666),(2,964890))]
    a,b = sites[0][1],sites[1][1]
    covariance = rational_covariance(adj(a)*a,adj(b)*b)
    levels,basis = model["eig"](covariance)
    assert levels[0]>0
    source = basis*acb_mat([[v.sqrt() if i==j else 0 for j in range(3)]
                            for i,v in enumerate(levels)])
    (up,U),(down,D) = [model["eig"](adj(y)*y) for y in (a*source*arb("7.06539940386823251"),b*source)]
    result = [v.sqrt() for v in up+down]+model["mixing_angles"](U,D)
    error = max((abs(x-y)/tol).upper() for x,y,tol in zip(result,model["ref"],model["tol"],strict=True))
    assert error < arb(".1232") < arb(1)/8
    print("固定有理二次多项式的六质量及四个 CKM 读数",[v.str(12) for v in result])
    print("距中心的最大参考半宽倍数",error.str(12))

    grams = [[adj(site)*site for site in sector] for sector in sites]
    for i,u in enumerate(grams[0],start=-1):
        for j,d in enumerate(grams[1],start=-1):
            if i==j==0:
                continue
            spectrum = rational_covariance(u,d).eig(algorithm="rump")
            assert all(lam.imag.contains(0) for lam in spectrum)
            assert sum(lam.real < arb("-3.88") for lam in spectrum)==2
            assert sum(lam.real > arb("4.46") for lam in spectrum)==1
            print("相邻编号偏移及协方差谱",i,j,[lam.real.str(9) for lam in spectrum])


def verify():
    orbit = run_path(str(Path(__file__).with_name("ch39-two-zero-orbit-condition.py")))
    model, up, down = orbit["calibration"]()
    certify_box(orbit,model,up,down)
    fixed_polynomial(model)


if __name__=="__main__":
    verify()

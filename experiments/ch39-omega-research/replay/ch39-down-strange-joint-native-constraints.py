"""任意右侧下型输运的双轻质量必要条件，及七项奇质量见证的距离禁区。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def adj(a):
    return a.conjugate().transpose()


def positive(a):
    for n in range(1,a.nrows()+1):
        assert acb_mat([[a[i,j] for j in range(n)] for i in range(n)]).det().real > 0


def verify(root):
    rank = run_path(str(root/"ch39-rank-two-generator-up-quark.py"))
    strange = run_path(str(root/"ch39-one-generator-gravity-up-strange.py"))
    train = run_path(str(root/"ch39-common-embedding-central-right-ckm.py"))["verify_nodes"]
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    spec = run_path(str(root/"ch39-adjacent-pvm-transport.py"))["spectrum"]
    n = rank["N"]+1
    _, labels, _ = train(root, ((n,".0294"),), ".1", (), ((6,7,10,13),(6,7,10,12)))
    ctx.prec = 384
    src, dest = [spec(shell,2,node)[1] for node in (n,n+1)]
    w = labels[n,2]
    c = acb_mat([[dest[i,j] for j in range(19)] for i in (14,8,7,6)])
    c4 = acb_mat([[c[i,j] for j in range(15,19)] for i in range(4)])
    coef = list(map(arb,(".7766",".98078",".97517",".7766")))
    positive(adj(c4)*c4 - acb_mat([[coef[i] if i==j else 0 for j in range(4)] for i in range(4)]))
    d,s = arb(".0000137"),arb(".000271")
    low = sorted(x.lower() for x in w)
    assert 1-d/min(coef)/low[10] > arb(".99972")
    assert 2-(d+s)/min(coef)/low[10] > arb("1.99435")
    omega = acb(-1,arb(3).sqrt())/2
    sqrt3 = arb(3).sqrt()
    chars = [[abs(c[1,j]+o*c[2,j]+o**2*c[3,j])**2/3 for j in range(19)]
             for o in (omega,omega.conjugate())]
    for first,other,mult,goal in ((0,1,"1.47",".00007359"),(1,0,"2.75296",".00009613")):
        t = arb(mult)
        order = sorted((max(arb(0),(t*a+b).lower()) for a,b in zip(chars[first],chars[other],strict=True)),
                       reverse=True)
        bound = sum(x*y for x,y in zip(low,order,strict=True)) - t*d
        assert bound > arb(goal)
        print("字符、另一个字符的条件下界", first+1, bound.lower().str(12))
    print("第十六至十九层系数", [v.str(7) for v in coef])
    print("已选两层进入最低十个源振幅的概率和下界",(2-(d+s)/min(coef)/low[10]).lower().str(12))

    x0 = strange["cayley"](root,rank["eye"])
    refine = run_path(str(root/"ch39-refine-two-gravity-images.py"))["high_images"]
    a,b,y,z = refine(root)
    predicted = x0*x0*x0
    err1, err2 = predicted*a-y,predicted*b-z
    mismatch = strange["frob"](acb_mat([[err1[i,0],err2[i,0]] for i in range(19)]))
    roots = x0.eig(algorithm="rump")
    dg = 4*mismatch
    assert min(abs(u**3-v**3).lower() for j,u in enumerate(roots) for v in roots[j+1:]) > 2*dg
    coefficient = min(abs(u*u+u*v+v*v).lower() for u in roots for v in roots)
    dx = 2*dg/coefficient
    assert (dg+3*dx**2+dx**3)/coefficient < dx and (6*dx+3*dx**2)/coefficient < 1
    assert dx < arb("2e-8")
    overlaps = adj(dest)*x0*src
    shifts = [arb(2)**(arb(k)/8) for k in range(-245,45)]
    inverses = [[sum(abs(overlaps[j,i])**2/(w[i]+g) for i in range(19)) for g in shifts] for j in range(19)]
    trivial = acb_mat([[1,0],[0,1/sqrt3],[0,1/sqrt3],[0,1/sqrt3]])
    for radius in (arb(".862"),arb("1.09")):
        rho = radius+dx
        overlap_sq = (1-rho*rho/2)**2
        assert rho < arb(2).sqrt()
        floors = [max([arb(0)]+[(overlap_sq/v-g).lower() for g,v in zip(shifts,row,strict=True)])
                  for row in inverses]
        cf = acb_mat([[c[i,j]*floors[j] for j in range(19)] for i in range(4)])
        mat = cf*adj(c)
        character_floors = [sum(f*q for f,q in zip(floors,p,strict=True)) for p in chars]
        tr = adj(trivial)*mat*trivial
        if radius < 1:
            positive(mat-rank["eye"](4)*d)
            print("排除轻下夸克的算子距离半径及两个字符下界", radius.str(6), [x.lower().str(12) for x in character_floors])
        else:
            positive(tr-rank["eye"](2)*s)
            print("联合入窗必须使用两条字符的算子距离半径",radius.str(7))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)

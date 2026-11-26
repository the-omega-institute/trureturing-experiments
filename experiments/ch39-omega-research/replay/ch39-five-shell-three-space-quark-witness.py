"""五个连续实际壳用三个复四站嵌入精确实现六个夸克汤川中心与完整 CKM。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


CENTRAL = 2_306_888_100_085


def weak_basis(a, b, coefficient, phase, diag, adj):
    vals, vec = (a+b*coefficient).eig(right=True, algorithm="rump")
    order = sorted(range(3), key=lambda i: vals[i].real.mid())
    assert all(vals[i].real < vals[j].real for i,j in zip(order, order[1:]))
    y = acb_mat([[vec[i,order[-1]]] for i in range(3)])
    y = y * (y[2,0].conjugate()/(abs(y[2,0])*(adj(y)*y)[0,0].real.sqrt()))
    residual = b*y - y*(adj(y)*b*y)[0,0].real
    z = residual * (-1/(adj(residual)*residual)[0,0].real.sqrt())
    x = acb_mat([[(y[(i+1)%3,0]*z[(i+2)%3,0]-y[(i+2)%3,0]*z[(i+1)%3,0]).conjugate()] for i in range(3)])
    basis = acb_mat([[vec[i,0] for vec in (x,y,z)] for i in range(3)])
    ratio = acb(1,phase)
    basis = basis*diag([ratio.conjugate()/ratio,1,1])
    assert (adj(basis)*basis).contains(diag([1]*3))
    return basis


def verify(root):
    fixed = run_path(str(root/"ch39-free-four-site-native-positive-spectrum.py"))
    two = run_path(str(root/"ch39-common-embedding-central-right-ckm.py"))
    solve = run_path(str(root/"ch39-complex-four-site-unit-quark-amplitudes.py"))
    ck = run_path(str(root/"ch39-shell-clock-ckm-blocks.py"))
    projections = run_path(str(root/"ch39-three-gram-source-flavor.py"))["projections"]
    ctx.prec = 320
    adj, diag = two["adj"], solve["diag"]
    mix = acb_mat(ck["standard_matrix"](*[arb(x).sin() for x in (".22704",".04275",".00377")], arb("1.139")))
    mass = [list(map(arb,row)) for row in ((".00000615",".00311",".8616"),(".0000135",".000268",".01401"))]
    a, b = [diag([v/row[-1] for v in row]) for row in mass]
    a = adj(mix)*a*mix
    u,v,intrinsic,rotation = two["common_embedding"](".00992624",".06181531",(".169285","-.022790",".747592","-.641817"))
    embedding = rotation*intrinsic
    assert (adj(embedding)*embedding).contains(diag([1]*3))
    params = ((0,".18779571","0",((7,10,13,16),(6,7,10,12))),
              (1,".30696321",".0126757",((6,7,10,13),(1,6,7,10))),
              (2,"-.03649316","-.0505816",((6,7,10,12),(6,7,10,12))))
    responses = {}
    for offset,t,phase,selections in params:
        number = CENTRAL+offset
        basis = weak_basis(a,b,arb(t),arb(phase),diag,adj)
        spectra = fixed["native"](root,number)
        projectors = []
        for qi,(charge,target,m,rows,idx) in enumerate(zip((-4,2),(a,b),mass,spectra,selections,strict=True)):
            expected = adj(basis)*target*basis*m[-1]
            four = two["extend_target"](expected,u,v)
            assert (adj(intrinsic)*four*intrinsic-expected).contains(acb_mat(3,3))
            _,weights = solve["reconstruct"](rows,four,rotation,[j-1 for j in idx])
            responses[number,charge] = weights
            physical = sum((w*row for w,row in zip(weights,rows,strict=True)),acb_mat(4,4))
            result = adj(embedding)*physical*embedding
            diff = result-expected
            assert diff.contains(acb_mat(3,3)) and max(abs(z).upper() for z in diff.entries())<arb("1e-40")
            eigen,proj = projections(result*result)
            assert all(e.overlaps(y*y) for e,y in zip(eigen,m,strict=True))
            projectors.append(proj)
            high = max(w.lower() for w in weights).union(max(w.upper() for w in weights))
            print("编号、荷、十九层峰值",number,charge,high.str(16))
        up,down = projectors
        assert all((up[i]*down[j]).trace().real.overlaps(abs(mix[i,j])**2) for i in range(3) for j in range(3))
        j = (up[0]*down[0]*up[1]*down[1]).trace().imag
        assert j>0 and j.overlaps((mix[0,0]*mix[1,1]*mix[0,1].conjugate()*mix[1,0].conjugate()).imag)
        print("编号、有符号 J",number,j.str(16))
    singles = ((-2,".1",(".00384",".0384",(".090494",".004868","-.054018",".994419")),((7,10,12,16),(7,10,12,16))),
               (-1,".74",(".0182",".113",("-.21668","-.53516",".71837","-.38807")),((1,10,16,18),(6,7,10,18))))
    frames = {CENTRAL+i: embedding for i in range(3)}
    for offset,t,geometry,picks in singles:
        n = CENTRAL+offset
        frame,weights,_ = two["verify_nodes"](root,((n,"0"),),t,geometry,picks)
        frames[n] = frame
        responses.update(weights)
    floor = arb("91.1876")/arb("1.22089e19")
    assert len(frames) == 5 and len(responses) == 10
    for number in range(CENTRAL-2,CENTRAL+3):
        for charge in (-4,2):
            weights = responses[number,charge]
            assert len(weights) == 19 and all(floor < w < 1 for w in weights)
    print("190 个夸克谱振幅均在共同截止下界与一之间", floor.str(16))
    return frames,responses


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)

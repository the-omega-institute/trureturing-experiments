"""有理矩阵见证：相同的原生谱线占据、四端口压缩与八站压缩。"""
from pathlib import Path
from runpy import run_path
from flint import acb_mat, arb, fmpq, fmpq_mat


def verify(root):
    graph = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    points = graph["geo"]["sites"]
    cycle = [points.index((y,z,x)) for x,y,z in points]
    run_path(str(root/"ch39-eleven-weak-eight-sites.py"))["check_cyclic_geometry"](graph,cycle)
    sites = run_path(str(root/"ch39-terminal-port-support.py"))["verify"](root)
    ports = (14,8,7,6)
    unit = fmpq_mat([[int(i==j) for j in range(19)] for i in range(19)])
    rotate = fmpq_mat([[int(cycle[j]==i) for j in range(19)] for i in range(19)])
    mean = (unit+rotate+rotate*rotate)/3
    flip = 2*mean-unit
    assert rotate**3==unit and mean**2==mean and flip**2==unit
    u = fmpq_mat([[int(i==14)] for i in range(19)])
    v = fmpq_mat([[fmpq(int(i in (0,2))-int(i in (3,5)),2)] for i in range(19)])
    assert mean*u==u and mean*v==fmpq_mat(19,1)
    plus,minus = ((u+s*v)*(u+s*v).transpose()/2 for s in (1,-1))
    common = fmpq_mat([[fmpq(int(i in (6,7,8) and j in (6,7,8)),3)
                         for j in range(19)] for i in range(19)])
    assert plus**2==plus and minus**2==minus and common**2==common
    assert plus*minus==common*plus==common*minus==fmpq_mat(19,19)
    assert flip*plus*flip==minus and flip*common*flip==common
    diff=plus-minus
    assert diff**2==plus+minus
    assert all(diff[i,i]==0 for i in range(19))
    for read in (sites,ports):
        assert all(diff[i,j]==0 for i in read for j in read)
    assert diff[14,0]==fmpq(1,2)
    print("四端口与八站压缩相同，极点与 0 站的归一化非对角元之差",diff[14,0])

    frames, _ = run_path(str(root/"ch39-cyclic-two-characters-nine.py"))["characters"](points)
    # 三条未调谐轨道在每个字符各贡献一个坐标。所有可调边的端点都在外部。
    quiet = {9,10,11,12,13,15,16,17,18}
    assert quiet.isdisjoint(set(graph["soft"]["star"]["orbit"])|{14})
    for i,j,d,_ in graph["geo"]["edges"]:
        assert d not in (724,804) or (i not in quiet and j not in quiet)
    minors = []
    for charge in (1,-4,2,-3,6):
        base = graph["laplacian"](charge,arb(1),phases={})
        for frame in frames:
            for k in range(frame.ncols()):
                support = {i for i in range(19) if abs(frame[i,k])>arb(".5")}
                assert bool(support & quiet) == (k in (3,4,5))
                if k in (3,4,5):
                    assert support <= quiet
            block = frame.conjugate().transpose()*base*frame
            bridge = acb_mat([block.tolist()[j][:3] for j in (3,4,5)])
            minors.append(abs(bridge.det()).lower())
    assert all(x.is_finite() and x>arb(".0001548") for x in minors)
    print("五荷十五个三阶连接子式的最小绝对值下界",min(minors).str(13))


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)

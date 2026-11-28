"""从两种实际弧相位的一阶谐波读取四坐标，在相邻两枚零点检验外推。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat

ROOT = Path(__file__).resolve().parent


def verify():
    phases = run_path(str(ROOT/"ch39-arc-relative-determinant.py"))["phases"]
    series = run_path(str(ROOT/"ch39-five-zero-split-rank.py"))["systems"]
    certify = run_path(str(ROOT/"ch39-two-split-global-roots.py"))["certify"]
    def features(code):
        arc = phases(acb.zeta_zero(code).imag)
        return [arb(1),arc[724].real,arc[724].imag,arc[804].real,arc[804].imag]
    train = tuple(range(2_345_664,2_345_669))
    held = (2_345_663,2_345_669)
    roots = [certify(system,seed,beta,logscale)
             for _,system,seed,beta,logscale in series(train)]
    matrix = arb_mat([features(n) for n in train])
    assert not matrix.det().contains(0)
    coefficients = matrix.solve(arb_mat(roots))
    for code,system,_,beta,logscale in series(held):
        predicted = (arb_mat([features(code)])*coefficients).entries()
        assert all(0 < u < 1 for u in predicted)
        readings,_,_,_,projection,target,_ = system
        residual = projection*readings(predicted)-target
        # 只用预测值作牛顿初值；以下独立校准用于量出预测的坐标误差。
        root = certify(system,[u.mid() for u in predicted],beta,logscale)
        roots.append(root)
        delta = [p-r for p,r in zip(predicted,root,strict=True)]
        peak = max(abs(u).lower() for u in delta).union(max(abs(u).upper() for u in delta))
        print("留出编号、四个坐标误差、最大误差、两项引力与两项规范交汇残差",
              code,[z.str(8) for z in delta],peak.str(9),[z.str(9) for z in residual.entries()],flush=True)
        lows = {held[0]:".00783",held[1]:".00616"}
        assert peak > arb(lows[code])
        interval = (".007831",".007833") if code == held[0] else (".006161",".006164")
        assert arb(interval[0]) < peak < arb(interval[1])
        assert abs(residual[3,0]) > arb(".032" if code == held[0] else ".019")
        assert all(not z.contains(0) for z in delta)
    # 七点重新联合校准的离散一致逼近。对偶单位球顶点至多支撑六点。
    all_features = [features(n) for n in (*train,*held)]
    distances = [[] for _ in range(4)]
    for support in combinations(range(7),6):
        row = [all_features[i] for i in support]
        cofactors = [(-1)**k*arb_mat(row[:k]+row[k+1:]).det() for k in range(6)]
        assert all(not h.contains(0) for h in cofactors)
        norm = sum(abs(h) for h in cofactors)
        for j,values in enumerate(distances):
            values.append(abs(sum(h*roots[i][j] for h,i in zip(cofactors,support,strict=True)))/norm)
    errors = [max(v.lower() for v in values).union(max(v.upper() for v in values))
              for values in distances]
    print("两弧一阶谐波在七枚零点自由重拟合的逐坐标最小最大误差",
          [v.str(9) for v in errors],flush=True)
    assert all(v > 0 for v in errors)
    assert arb(".00034922") < errors[3] < arb(".00034923")


if __name__ == "__main__":
    verify()

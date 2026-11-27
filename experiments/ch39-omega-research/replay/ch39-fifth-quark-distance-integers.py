"""以整数惯性与有理柯西界认证第五壳轻夸克窗所需的第四步位移。"""
from fractions import Fraction
from pathlib import Path
from runpy import run_path


def verify(root):
    rows = [line.split() for line in (root/"ch39-fifth-quark-distance-cells.tsv").read_text().splitlines()]
    scale = 10**14
    sites = (14,8,7,6)
    assert rows[:3]==[["scale",str(scale)],["shell","2306888100090"],["sites",*map(str,sites)]]
    assert len(rows)==51
    inertia = run_path(str(root/"ch39-fourth-shell-strange-integers.py"))
    shifts = [Fraction(2)**k*m/4 for k in range(-32,10) for m in (4,5,6,7)]
    tests = ((-4,2,("0.0344","0.00000629"),("0.0474","0.00316")),
             (2,1,("0.0331","0.0000137"),("0.0261","0.000271")))
    for block, (charge,power,*windows) in enumerate(tests):
        data = rows[3+24*block:3+24*(block+1)]
        assert data[0][:2]==["train",str(charge)] and len(data[0])==21
        source = [Fraction(int(v),scale) for v in data[0][2:]]
        assert all(0<=v<1 for v in source)
        resolvents = []
        for i,row in enumerate(data[1:20]):
            assert row[:3]==["prob",str(charge),str(i)] and len(row)==22
            probabilities = [Fraction(int(v),scale) for v in row[3:]]
            assert min(probabilities)>=0 and sum(probabilities)>=1
            resolvents.append([sum(p/(t+s) for p,t in zip(probabilities,source,strict=True)) for s in shifts])
        basis=[]
        for site,row in zip(sites,data[20:],strict=True):
            assert row[:3]==["basis",str(charge),str(site)] and len(row)==41
            values=list(map(int,row[3:]))
            basis.append(list(zip(values[::2],values[1::2],strict=True)))
        for j,(radius,ceiling) in enumerate(windows):
            distance = power*(Fraction(radius)+2*Fraction("1.28e-10"))
            overlap = 1-distance**2/2
            assert overlap>0
            weights = [max(0,max(int((overlap**2/r-s)*scale) for r,s in zip(row,shifts,strict=True)))
                       for row in resolvents]
            assert max(weights)<=scale
            cells={}
            for i in range(4):
                for k in range(i,4):
                    pairs=list(zip(basis[i],basis[k],weights,strict=True))
                    real=sum(w*(a*c+b*d) for (a,b),(c,d),w in pairs)
                    imag=sum(w*(b*c-a*d) for (a,b),(c,d),w in pairs)
                    cells[i,k]=(real,imag)
            # 谱框格点误差的算子范数 < 13/scale，压缩矩阵误差 < 27/scale。
            shift = (Fraction(ceiling)+Fraction(27,scale))*scale**3
            assert shift.denominator==1
            count=inertia["negative_count"](inertia["realify"](4,cells),int(shift))
            assert count<=2*j
            print("荷、三代能序、排除半径、质量窗上端、实化负惯性",charge,j+1,radius,ceiling,count)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)

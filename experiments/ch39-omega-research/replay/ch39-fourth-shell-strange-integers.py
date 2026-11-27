"""核验第四壳奇夸克下界的有限整数证书。仅使用 Python 标准库。"""
from pathlib import Path


SCALE = 10**12
CUT_SCALE = 10**9
CUTS = {
    7: (197,1373,7057,97402,5905026,20410143,69023069),
    12: (10059,10436,165402,49180580,97172003,118011563,147538478,153368430,
         153770386,153773264,153774798,273164596),
}


def require(condition, detail):
    if not condition:
        raise ValueError(detail)


def read(path):
    rows = [line.split("\t") for line in path.read_text(encoding="utf-8").splitlines()]
    header = [["scale",str(SCALE)],["entry_error","1"],["shells","2306888100086","2306888100089"],
              ["sites","14","8","7","6"],["trial","1","0","-2","5"]]
    require(rows[:5]==header,"格距、壳、站号或试探矢不符")
    groups = []
    for row,tag,n in zip(rows[5:7],("trivial","complement"),(7,12),strict=True):
        ids = list(map(int,row[1:]))
        require(row[0]==tag and len(ids)==n and sorted(set(ids))==ids,"块内谱线表无效")
        groups.append(set(ids))
    require(groups[0]=={1,6,7,10,12,16,19} and groups[0]|groups[1]==set(range(1,20)),
            "两块未给出指定的十九层分割")
    blocks = {n:{} for n in CUTS}
    rays = {j:{} for j in range(3)}
    for tag,*data in rows[7:]:
        if tag=="cell":
            n,i,j,a,b = map(int,data)
            require(n in blocks and 0<=i<=j<n and (i,j) not in blocks[n] and (i!=j or b==0),
                    "格点坐标重复或无效")
            blocks[n][i,j] = a,b
        elif tag=="ray":
            q,j,v = map(int,data)
            require(q in rays and j not in rays[q] and 0<=v<=SCALE,"射线整数下界无效")
            rays[q][j] = v
        else:
            raise ValueError("未知行类")
    require(all(len(blocks[n])==n*(n+1)//2 for n in CUTS),"上三角不完整")
    require(set(rays[0])==groups[0] and set(rays[1])|set(rays[2])==groups[1] and
            len(rays[1])==len(rays[2])==6,"射线与两个源块不匹配")
    return blocks,rays


def realify(n,cells):
    matrix = [[0]*(2*n) for _ in range(2*n)]
    for (i,j),(a,b) in cells.items():
        for k,l,z in ((i,j,b),(j,i,-b)):
            matrix[k][l] = matrix[k+n][l+n] = a
            matrix[k+n][l],matrix[k][l+n] = z,-z
    return matrix


def negative_count(matrix,shift):
    # 第 j 个无换行 Bareiss 枢轴为 j 阶前导主子式。相邻主子式异号即一个负惯性。
    a = [row.copy() for row in matrix]
    for i in range(len(a)):
        a[i][i] -= shift
    previous = 1
    negative = 0
    for k in range(len(a)):
        pivot = a[k][k]
        require(pivot!=0,"给定分隔处有零前导主子式")
        negative += (pivot<0)!=(previous<0)
        for i in range(k+1,len(a)):
            for j in range(i,len(a)):
                value = pivot*a[i][j]-a[i][k]*a[k][j]
                quotient,remainder = divmod(value,previous)
                require(remainder==0,"Bareiss 除法未整除")
                a[i][j]=a[j][i]=quotient
        previous = pivot
    return negative


def verify(path):
    blocks,rays = read(path)
    for n,cuts in CUTS.items():
        m = realify(n,blocks[n])
        counts = [negative_count(m,cut*(SCALE//CUT_SCALE)+2*n) for cut in cuts]
        require(counts==list(range(0,2*n,2)),f"{n} 维源谱的分隔惯性错误")
        print("块维数与移位后的负惯性",n,counts)
    denominator = SCALE*CUT_SCALE
    lower = []
    for sector,row in rays.items():
        weights = CUTS[7 if sector==0 else 12]
        coeff = sorted(row.values(),reverse=True)
        coeff.extend([0]*(len(weights)-len(coeff)))
        lower.append(sum(a*b for a,b in zip(weights,coeff,strict=True)))
    for sector,value in enumerate(lower):
        print("射线组与逆序配对下界的分子、分母",sector,value,denominator)
    require(min(lower)*10**8>27621*denominator,"未能核验奇质量 > 0.00027621")
    print("奇夸克 Yukawa 严格大于 27621/100000000")
    return lower,denominator


if __name__=="__main__":
    verify(Path(__file__).with_name("ch39-fourth-shell-strange-cells.tsv"))

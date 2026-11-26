"""用原 96 位有理 Gram 格计算所有主 Schur 残向量对的交叉谱负量。"""
from fractions import Fraction as Q
from itertools import combinations, permutations
from math import isqrt
from pathlib import Path
from runpy import run_path
import sys


ZERO = (Q(0),Q(0))
ONE = (Q(1),Q(1))


def add(a,b):
    return a[0]+b[0],a[1]+b[1]


def neg(a):
    return -a[1],-a[0]


def mul(a,b):
    ends = (a[0]*b[0],a[0]*b[1],a[1]*b[0],a[1]*b[1])
    return min(ends),max(ends)


def det(matrix):
    n = len(matrix)
    result = ZERO
    for p in permutations(range(n)):
        product = ONE
        for i,j in enumerate(p):
            product = mul(product,matrix[i][j])
        inverse = sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
        result = add(result,neg(product) if inverse%2 else product)
    return result


def coefficients(g,i,s):
    # 乘正数 det G_ss 后，主残向量的系数只需乘加和精确有理端点。
    m = [[g[a][b] for b in s] for a in s]
    c = [ZERO]*4
    c[i] = det(m)
    assert c[i][0]>0
    for k,j in enumerate(s):
        replaced = [[g[a][i] if col==k else value for col,value in enumerate(row)]
                    for a,row in zip(s,m)]
        c[j] = neg(det(replaced))
    return c


def dot(g,a,b):
    value = ZERO
    for i in range(4):
        for j in range(4):
            value = add(value,mul(mul(a[i],g[i][j]),b[j]))
    return value


def sqrt_bounds(bounds):
    scale = 10**14
    a,b = bounds
    low = isqrt(a.numerator*scale**2//a.denominator)
    high = isqrt(b.numerator*scale**2//b.denominator)+1
    return Q(low,scale),Q(high,scale)


def verify(details=False):
    source = run_path(str(Path(__file__).with_name("ch39-finite-rh-pivot-mass.py")))
    g = [[source["read"](i,j) for j in range(4)] for i in range(4)]
    assert det(g)[0]>0
    cols = {(i,s): coefficients(g,i,s) for i in range(4) for n in range(4)
            for s in combinations([j for j in range(4) if j!=i],n)}
    norms = {key:dot(g,c,c) for key,c in cols.items()}
    assert min(value[0] for value in norms.values())>0
    groups = {name:[] for name in ("原向量","同一非空消元集","原向量与残向量","不同非空消元集")}
    zeros = 0
    for a,b in combinations(cols,2):
        i,s = a
        j,t = b
        if {i,*s}<=set(t) or {j,*t}<=set(s):
            zeros += 1
            if details:
                print(a,b,"严格正交，分母为零")
            continue
        scalar = dot(g,cols[a],cols[b])
        assert scalar[0]>0 or scalar[1]<0,(a,b)
        denominator = mul(scalar,scalar)
        numerator = mul(norms[a],norms[b])
        tv = sqrt_bounds((numerator[0]/denominator[1],numerator[1]/denominator[0]))
        low,high = [(e-1)/2 for e in tv]
        if not s and not t:
            name = "原向量"
        elif s==t:
            name = "同一非空消元集"
        elif not s or not t:
            name = "原向量与残向量"
        else:
            name = "不同非空消元集"
        groups[name].append((low,high,a,b))
        if details:
            print(a,b,"最大负量的有理包围",low,high)
    assert zeros==108
    # 左右比较都避开物理证书舍入端点。
    costs = [(Q("1.86547379546839"),Q("1.86547379546841")),
             (Q("43.42807"), Q("43.42808")),
             (Q("1308843.19304049"),Q("1308843.19305"))]
    results = [(6,0,0,0),(18,0,0,0),(64,64,28,0),(300,180,32,0)]
    for (name,rows),expected in zip(groups.items(),results,strict=True):
        passed = []
        for left,right in costs:
            assert all(high<left or low>right for low,high,_,_ in rows)
            passed.append(sum(low>right for low,_,_,_ in rows))
        assert tuple([len(rows),*passed])==expected
        maximum = max(rows,key=lambda row:row[0])
        assert all(row==maximum or row[1]<maximum[0] for row in rows)
        lo,hi,a,b = maximum
        print(name,len(rows),"最大负量",f"{float(lo):.9f}",f"{float(hi):.9f}",a,b,
              "达到三项预算的对数",*passed)
    largest = max(q[1] for rows in groups.values() for q in rows)
    assert Q("3844.09") < largest < Q("3844.11")
    print("相互正交",zeros,"全部可归一对",sum(len(v) for v in groups.values()))
    return groups


if __name__=="__main__":
    verify("--pairs" in sys.argv[1:])

"""三编号共同味多项式的标准库整数秩证书。"""
from fractions import Fraction as F
from itertools import combinations, permutations, product
from math import isqrt
from pathlib import Path
from runpy import run_path

BITS = 130
WORK = 224
Q = 1 << WORK
RECORD = [
    (1352417259997240471391108987852029198877, -153756244291461525726507447520669293374, 1360782735945456175244706985904689090639, 30720894353160048243762503284267050677, 3350256319),
    (1361054906855856166645865411976941798214, -14246695078361899061081456258327747408, 1347965787568008940232591183047393878419, 188842959474841214338592181530113724041, 3350256319),
    (1360851203709981296798761249366158650355, 27521430888565978790423125755695855174, 1340571242256911017295262399953226071871, 235673868366481224593889231938708609539, 3350256319),
]


def dot(x,y):
    return sum(a*b for a,b in zip(x,y,strict=True))


def geometry():
    outer=sorted(set(permutations((1,5,29))))
    inner=sorted(set(permutations((13,13,23))))
    sites=outer+inner+sorted(x for x in product(range(30),repeat=3)
                             if dot(x,x)==867 and x not in outer+inner)
    pole=(17,17,17)
    positions=[sites.index(x) for x in (pole,(23,13,13),(13,23,13),(13,13,23))]
    assert len(sites)==19 and positions==[14,8,7,6]
    edges=[]
    for i,j in combinations(range(19),2):
        x,y=sites[i],sites[j]
        dist=2*(867-dot(x,y))
        if dist>804:
            continue
        det=sum(pole[k]*(x[(k+1)%3]*y[(k+2)%3]-x[(k+2)%3]*y[(k+1)%3]) for k in range(3))
        a=867*(867+dot(pole,x)+dot(x,y)+dot(y,pole))
        b=-det
        denom=a*a+867*b*b
        re=isqrt(Q*Q*a*a//denom)*(1 if a>=0 else -1)
        im=isqrt(Q*Q*867*b*b//denom)*(1 if b>=0 else -1)
        if dist==200:
            re,im=-re,-im
        edges.append((i,j,dist,(re,im)))
    assert len(edges)==135
    return positions,edges


def mul(z,w):
    a,b=z
    c,d=w
    return (a*c-b*d)//Q,(a*d+b*c)//Q


def gauss_matrix(record,charge,edges):
    re=[[0]*19 for _ in range(19)]
    im=[[0]*19 for _ in range(19)]
    u,v,s,t,delta=[(2*k+1)<<(WORK-BITS-1) for k in record]
    for i,j,d,geom in edges:
        weight=Q-delta if 14 in (i,j) and (i in (0,3,4) or j in (0,3,4)) else Q
        unit=(Q,0)
        if 14 not in (i,j):
            clock={724:(u,v),804:(s,t)}.get(d,(Q,0))
            single=mul(geom,clock)
            if charge<0:
                single=(single[0],-single[1])
            for _ in range(abs(charge)):
                unit=mul(unit,single)
        h,k=mul((-weight,0),unit)
        re[i][j]=re[j][i]=h
        im[i][j],im[j][i]=k,-k
        re[i][i]+=weight
        re[j][j]+=weight
    for i in range(19):
        re[i][i]-=17*Q
    return [[(re[i][j]//16,im[i][j]//16) for j in range(19)] for i in range(19)]


def moment_rows(matrix,positions):
    y=[[(Q*int(i==j),0) for j in positions] for i in range(19)]
    cols=[]
    for degree in range(36):
        p,r=positions[:2]
        cols.append([y[p][0][0],y[r][1][0],*y[p][1],*y[r][2]])
        if degree==35:
            break
        updated=[]
        for row in matrix:
            vals=[]
            for j in range(4):
                real=imag=0
                for (a,b),v in zip(row,y,strict=True):
                    c,d=v[j]
                    real+=a*c-b*d
                    imag+=a*d+b*c
                vals.append((real//Q,imag//Q))
            updated.append(vals)
        y=updated
    return [list(row) for row in zip(*cols,strict=True)]


def inverse_norm(matrix):
    n=len(matrix)
    augmented=[row+[int(i==j) for j in range(n)] for i,row in enumerate(matrix)]
    previous=1
    for i in range(n):
        p=next(k for k in range(i,n) if augmented[k][i])
        augmented[i],augmented[p]=augmented[p],augmented[i]
        row=augmented[i]
        pivot=row[i]
        for k in range(n):
            if k==i:
                continue
            r=augmented[k]
            left=r[i]
            for j in range(2*n):
                if j!=i:
                    quotient,remainder=divmod(pivot*r[j]-left*row[j],previous)
                    assert remainder==0
                    r[j]=quotient
            r[i]=0
        previous=pivot
    signed=augmented[-1][-n-1]
    assert signed!=0
    invers=[r[n:] for r in augmented]
    for i,row in enumerate(invers):
        for j,col in enumerate(zip(*matrix,strict=True)):
            assert dot(row,col)==signed*int(i==j)
    return F(Q*max(sum(abs(v) for v in row) for row in invers),abs(signed))


def gram_injection(root):
    # 法向量首分量取一。n1<0<n2 使其同时具有平凡分量、横向分量。
    lib=run_path(str(root/"ch39-finite-rh-pivot-mass.py"))
    read=lib["read"]
    x0,x1=read(1,0)
    y0,y1=read(2,0)
    z0,_=read(2,1)
    assert 0<x0<x1<1 and 0<y0<y1
    # n2 的符号等于 x*z-y；三个源独立由其 Gram 主子式证实。
    assert x0*z0>y1
    assert all(lib["determinant"](range(k))[0]>0 for k in (2,3,4))


def verify(root):
    positions,edges=geometry()
    matrix=[]
    for record in RECORD:
        for charge in (-4,2):
            matrix.extend(moment_rows(gauss_matrix(record,charge,edges),positions))
    norm=inverse_norm(matrix)
    assert norm<2**108
    # 逐边误差、标度壳行和和四列传播余数均为两整数之比。
    epsilon=F(1,2**BITS)
    h=F(1,Q)
    assert 4*(1+epsilon+5*h)**3<5
    step=(18*epsilon+18*(5*epsilon+37*h))/16+38*h
    assert step<10*epsilon+128*h
    assert F(19,16)+step<F(5,4)
    growth=sum(k*F(5,4)**(k-1) for k in range(1,36))
    def defect(bits):
        return norm*growth*(10*F(1,2**bits)+F(136,Q))
    first=next(k for k in range(20,BITS+1) if defect(k)<1)
    print("逆矩阵行和小于 2^108，中心半径误差界通过的首个深度",first)
    def ceil_bound(value,scale):
        return (value.numerator*scale+value.denominator-1)//value.denominator
    print("130 位格的预条件缺陷上界（百万分之一）",ceil_bound(defect(BITS),10**6))
    assert first==130 and defect(BITS)<F(520523,10**6)
    gram_injection(root)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)

"""用纯整数行列式穷举七条平凡射线，并检验两个非平凡字符。"""
from itertools import permutations
from pathlib import Path


SCALE = 10**16
LOWER = 71789*SCALE**2//10**10
UPPER = 71792*SCALE**2//10**10
RAYS = (1,4,7,10,12,16,18)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read(path):
    rows = [line.split("\t") for line in path.read_text(encoding="utf-8").splitlines()]
    require(rows[:3]==[["scale",str(SCALE)],["shells","2306888100087","2306888100088"],
                       ["sites","14","8","7","6"]], "外包表的缩放、壳或站位错误")
    weights, rays, chars = {}, {}, {1:{},2:{}}
    for tag,*fields in rows[3:]:
        values = list(map(int,fields))
        if tag=="weight":
            j,lo,hi = values
            require(j not in weights and 0<=lo<=hi<=SCALE, "源权行重复或外包无效")
            weights[j] = (lo,hi)
        elif tag=="ray":
            j,*ends = values
            require(j not in rays and len(ends)==8, "平凡射线重复或列数错误")
            ranges = list(zip(ends[::2],ends[1::2],strict=True))
            require(all(-SCALE<=l<=h<=SCALE for l,h in ranges), "射线外包端点无效")
            rays[j] = ranges
        elif tag=="character":
            k,j,lo,hi = values
            require(k in chars and j not in chars[k] and 0<=lo<=hi<=SCALE, "字符行无效")
            chars[k][j] = (lo,hi)
        else:
            raise ValueError("未知的外包表行类")
    require(set(weights)==set(range(1,20)) and tuple(sorted(rays))==RAYS, "源层或平凡层编号不完整")
    levels = [*rays,*chars[1],*chars[2]]
    require(len(chars[1])==len(chars[2])==6 and sorted(levels)==list(range(1,20)), "三个字符的层未分割目的谱")
    return weights,[rays[i] for i in RAYS],chars


def contraction(weights, rays):
    """矩阵两项对角、非对角实部和虚部的整数分子区间，分母为 SCALE²。"""
    return [(sum(w*lo for w,(lo,_) in zip(weights,col,strict=True)),
             sum(w*hi for w,(_,hi) in zip(weights,col,strict=True))) for col in zip(*rays,strict=True)]


def square_low(interval):
    a,b = interval
    return 0 if a<=0<=b else min(a*a,b*b)


def positive_above(matrix, level):
    (a,_),(d,_),re,im = matrix
    a -= level
    d -= level
    determinant = a*d-max(x*x for x in re)-max(x*x for x in im)
    return a>0 and d>0 and determinant>0


def negative_direction(matrix, level):
    (a,b),(c,d),re,im = matrix
    a,b,c,d = a-level,b-level,c-level,d-level
    det_upper = max(a*c,a*d,b*c,b*d)-square_low(re)-square_low(im)
    return b<0 or d<0 or det_upper<0


def verify(path):
    weights,rays,chars = read(path)
    low = sorted(v[0] for v in weights.values())
    high = sorted(v[1] for v in weights.values())
    for j,p in chars.items():
        coeff = sorted((v[0] for v in p.values()), reverse=True)
        value = sum(w*z for w,z in zip(low[:6],coeff,strict=True))
        require(value>LOWER, f"字符 {j} 没有超过目标下界")
        print("字符、逆序配对积的整数下界",j,value)
    count = 0
    upper_witness = None
    for perm in permutations(range(7)):
        mat_low = contraction([low[j] for j in perm],rays)
        require(positive_above(mat_low,LOWER), f"排列 {perm} 的下界核验失败")
        if upper_witness is None:
            mat_high = contraction([high[j] for j in perm],rays)
            if negative_direction(mat_high,UPPER):
                upper_witness = tuple(j+1 for j in perm)
        count += 1
    require(count==5040 and upper_witness is not None, "未完成穷举或未找到上界排列")
    print("严格通过下界的排列数",count)
    print("可达上界的源权升序秩，按目的层排列",upper_witness)
    print("最轻上型最优值的严格下、上界", "71789/10000000000", "71792/10000000000")
    return upper_witness


if __name__=="__main__":
    verify(Path(__file__).with_name("ch39-next-shell-seven-ray-bounds.tsv"))

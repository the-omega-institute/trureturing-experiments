"""Find the shortest norm-five quaternion words entering the low-energy coupling box."""
from fractions import Fraction as Q
from math import isqrt

CENTERS = tuple(map(Q, (".467453",".63811","1.0583")))
BANDS = [(c*Q(999,1000),c*Q(1001,1000)) for c in CENTERS]
GENERATORS = [(1,2,0,0),(1,-2,0,0),(1,0,2,0),(1,0,-2,0),(1,0,0,2),(1,0,0,-2)]


def multiply(p,q):
    a,b,c,d = p
    w,x,y,z = q
    return (a*w-b*x-c*y-d*z, a*x+b*w+c*z-d*y,
            a*y+c*w+d*x-b*z, a*z+d*w+b*y-c*x)


def factor(q, length):
    letters = []
    for _ in range(length):
        choices = []
        for j in range(6):
            product = multiply(q,GENERATORS[j^1])
            if all(v%5 == 0 for v in product):
                choices.append((j,tuple(v//5 for v in product)))
        assert len(choices) == 1
        j,q = choices[0]
        letters.append(j)
    assert q == (1,0,0,0)
    letters.reverse()
    assert all((a^1) != b for a,b in zip(letters,letters[1:]))
    return letters


def ceil_rational(x):
    return -((-x.numerator)//x.denominator)


def hits(depth, bands=BANDS):
    norm = 5**depth
    low = Q(norm)/(1+sum(b*b for a,b in bands))
    high = Q(norm)/(1+sum(a*a for a,b in bands))
    start = isqrt(low.numerator//low.denominator)
    stop = isqrt(high.numerator//high.denominator)
    result = []
    for a in range(max(1,start),stop+1):
        if a%2 == 0 or a*a < low:
            continue
        endpoints = [(ceil_rational(a*l),int(a*r)) for l,r in bands]
        (b0,b1),(c0,c1),(d0,d1) = endpoints
        for b in range(b0+b0%2,b1+1,2):
            for c in range(c0+c0%2,c1+1,2):
                remaining = norm-a*a-b*b-c*c
                if not d0*d0 <= remaining <= d1*d1:
                    continue
                d = isqrt(remaining)
                if d*d == remaining and d%2 == 0 and any(v%5 for v in (a,b,c,d)):
                    q = (a if a%4 == 1 else -a,b,c,d)
                    letters = factor(q,depth)
                    back = (1,0,0,0)
                    for j in letters:
                        back = multiply(back,GENERATORS[j])
                    assert back == q
                    result.append((q,letters))
    return result


if __name__ == "__main__":
    assert all(not hits(d) for d in range(16))
    result = hits(16)
    assert len(result) == 4
    print("Shortest length 16, absolute coupling triples",len(result))
    for quat,letters in result:
        print(quat, "word",",".join(str(j) for j in letters))

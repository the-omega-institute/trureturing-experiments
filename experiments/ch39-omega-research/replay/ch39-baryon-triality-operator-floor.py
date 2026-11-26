"""Grassmann monomials witnessing the leading triple-baryon operators."""
from itertools import permutations, product


def sign(perm):
    return (-1)**sum(perm[i] > perm[j] for i in range(len(perm)) for j in range(i+1,len(perm)))


def monomial_coefficient(word, target):
    if set(word) != set(target):
        return 0
    return sign([target.index(x) for x in word])


def baryon(kind, target):
    result = [0,0]
    for colors in permutations(range(3)):
        for spin in range(2):
            for free in range(2):
                word = [("u",colors[0],spin), ("d",colors[1],1-spin), (kind,colors[2],free)]
                result[free] += sign(colors)*(-1)**spin*monomial_coefficient(word,target)
    return result


def qqql():
    target = [("Q",0,0,0),("Q",1,1,1),("Q",2,0,0),("L",0,1,1)]
    total = 0
    for colors in permutations(range(3)):
        for a,c,s,r in product(range(2),repeat=4):
            word = [("Q",colors[0],a,s),("Q",colors[1],1-a,1-s),
                    ("Q",colors[2],c,r),("L",0,1-c,1-r)]
            total += sign(colors)*(-1)**(a+c+s+r)*monomial_coefficient(word,target)
    return total


def verify():
    n1 = baryon("d",[("u",0,0),("d",1,1),("d",2,0)])
    n2 = baryon("d",[("u",0,1),("d",1,0),("d",2,1)])
    p3 = baryon("u",[("u",0,0),("d",1,1),("u",2,0)])
    assert n1[0]*n2[1]*p3[0] and n1[1] == n2[0] == p3[1] == 0
    # Six times hypercharge for 4 u^c, 5 d^c, e^c. Shell charge uses signed lifts.
    assert 4*(-4)+5*2+6 == 0 and (4*(-1)+5+(-1))%3 == 0
    assert qqql() != 0 and ((-1)*3)%3 == 0
    print("d=15 monomial coefficient", n1[0]*n2[1]*p3[0])
    print("d=18 three-family monomial coefficient", qqql()**3)


if __name__ == "__main__":
    verify()

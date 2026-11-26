"""二进制加法链的规范单态积分与五态费米子行列式。"""
from flint import arb, ctx
from sympy import Matrix, diag, eye, factor, factorial, symbols


def addition_chain(exponent):
    powers = [1]
    parents = []
    for bit in bin(exponent)[3:]:
        previous = len(powers)-1
        parents.append((previous, previous))
        powers.append(2*powers[previous])
        if bit == "1":
            parents.append((len(powers)-1, 0))
            powers.append(powers[-1]+1)
    weights = [0]*len(powers)
    weights[-1] = 1
    for child in range(len(parents), 0, -1):
        for parent in parents[child-1]:
            weights[parent] += weights[child]
    return powers, weights


def verify():
    powers, weights = addition_chain(23226)
    assert len(powers) == 23 and powers[-1] == weights[0] == 23226
    square_sum = sum(k*k for k in weights[1:])
    assert square_sum == 110729881 and max(weights[1:]) == 8192
    # 五个对角质量来自四条 B0 和一处裸质量，五条循环边中一条携带单态。
    z, t = symbols("z t", nonzero=True)
    cycle = Matrix([[0,0,0,0,z], [1,0,0,0,0], [0,1,0,0,0],
                    [0,0,1,0,0], [0,0,0,1,0]])
    mass = diag(2,3,5,7,11) + diag(1,4,6,8,9)*cycle
    gram = mass.T.subs(z, 1/z)*mass
    determinant = (t*eye(5)+gram).det(method="berkowitz")
    product = 2*3*5*7*11*1*4*6*8*9
    assert factor(determinant-determinant.subs(z,1)-product*(z+1/z-2)) == 0
    assert -16+6 == -10 and -16+4 == -12

    ctx.prec = 150
    kappa = arb(10)**9
    first_cos = kappa.bessel_i(1, scaled=True)/kappa.bessel_i(0, scaled=True)
    lower = 1-square_sum*(1-first_cos)
    assert lower > arb(".94463")
    # 在 kappa=1，仅最大的一次折叠已足以给出此上界。
    order = max(weights[1:])
    upper = (arb(1)/2)**order/arb(int(factorial(order)))
    assert upper < arb("1e-30969")
    print("规范单态数",len(powers)-1,"残差权平方和",square_sum,
          "强刚度一阶谐波传递下界",lower.str(14),"单位刚度上界",upper.str(9))


if __name__ == "__main__":
    verify()

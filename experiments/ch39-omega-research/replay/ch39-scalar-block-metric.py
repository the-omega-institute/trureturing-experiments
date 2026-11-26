"""Affine scalar Dirichlet data on the two signed Regge witnesses."""
from collections import defaultdict
from itertools import combinations, permutations

from flint import arb, arb_mat, ctx, fmpq, fmpq_mat

ctx.prec = 160
orders = tuple(permutations(range(4)))
zero = lambda n,m: fmpq_mat(n,m)


def simplex(order, parity):
    nodes = [[int(j in order[:i]) for j in range(4)] for i in range(5)]
    P = fmpq_mat(nodes[1:])
    q = {(i,i): fmpq(0) for i in range(5)}
    for i,j in combinations(range(5),2):
        d = j-i
        v = fmpq(d*(5-d),5)
        if parity:
            v *= (1+fmpq(parity*(-1)**i,24))*(1+fmpq(parity*(-1)**j,24))
        elif d == 1 and order[i] == 0:
            v += fmpq(1,3)
        q[i,j] = q[j,i] = v
    gram = fmpq_mat([[(q[0,i]+q[0,j]-q[i,j])/2 for i in range(1,5)] for j in range(1,5)])
    grad = fmpq_mat([[-1]*4]+[[int(i==j) for j in range(4)] for i in range(4)])
    local = grad*gram.inv()*grad.transpose()
    return gram.det(), local*fmpq_mat(nodes), P.transpose()*gram.inv()*P


def verify():
    for parities in ((-1,1),(0,)):
        tensors = defaultdict(lambda:zero(4,4))
        residual = defaultdict(lambda:zero(1,4))
        volume = arb(0)
        for parity in parities:
            for order in orders:
                det,force,inverse = simplex(order,parity)
                assert det > 0
                tensors[det] += inverse
                volume += arb(det).sqrt()/(24*len(parities))
                for i in range(5):
                    key = (parity*(-1)**i,det)
                    residual[key] += fmpq_mat([[force[i,j] for j in range(4)]])
        assert all(row == zero(1,4) for row in residual.values())
        D = arb_mat(4,4)
        for det,total in tensors.items():
            D += arb_mat(total)*(arb(det).sqrt()/(24*len(parities)))
        assert D.det() > 0
        ratio = D.det().sqrt()/volume
        lo,hi = (100558,100559) if len(parities)==2 else (101958,101959)
        assert fmpq(lo,100000) < ratio and ratio < fmpq(hi,100000)
        print(parities,"Dirichlet metric volume / microscopic volume",ratio)


if __name__ == "__main__":
    verify()

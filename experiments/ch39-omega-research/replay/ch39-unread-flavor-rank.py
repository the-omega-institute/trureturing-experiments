"""Five independent operators in both the continuum and rational shell flavor implementations."""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, nmod_mat

root = Path(__file__).parent
analytic = run_path(str(root/"ch39-fixed-local-flavor.py"))
finite = run_path(str(root/"ch39-rational-local-flavor.py"))
indices = (10815, 10965, 11177, 11431, 11831)
pair_indices = (11943,12019,11848,10898,11734,10965,11831,11468,11962,11905,
                10993,11639,11829,11620,11601)


def pair_coordinates(block):
    a,b,c,d = block
    return [1,a,b,c,d,a*d-b*c,a*a,a*b,b*b,2*a*c,a*d+b*c,2*b*d,c*c,c*d,d*d]


def verify():
    basis = analytic["basis"]
    varying = [arb(n)/2048+arb(0, arb(2)**-28) for n in analytic["codes"]]
    rows = []
    for n in indices:
        response = analytic["native"]["plane_response"](acb.zeta_zero(n).imag)
        w = analytic["mixing"](response, varying)
        block = basis.transpose()*w*basis
        rows.append([1, *block.entries()])
    determinant = abs(acb_mat(rows).det())
    assert determinant > arb("8e-10")

    prime = 65537
    def reduce_rational(q):
        return int(q.p)*pow(int(q.q) % prime, -1, prime) % prime
    finite_rows = []
    table = dict(finite["records"]())
    for n in indices:
        w = finite["flavor_operation"](table[n])
        finite_rows.append([1]+[(reduce_rational(w[i,j])+256*reduce_rational(w[i+3,j])) % prime
                                for i in range(2) for j in range(2)])
    residue = int(nmod_mat(finite_rows, prime).det())
    assert (256*256+1) % prime == 0 and residue == 21468
    print("Analytic five-operator determinant modulus", determinant.str(5))
    print("Gaussian-rational minor mod 65537:", residue, "exact Kraus rank 5")

    center = [arb(n)/2048 for n in analytic["codes"]]
    analytic_pairs, finite_pairs = [], []
    for n in pair_indices:
        response = analytic["native"]["plane_response"](acb.zeta_zero(n).imag)
        w = analytic["mixing"](response, center)
        block = basis.transpose()*w*basis
        analytic_pairs.append(pair_coordinates(block.entries()))
        v = finite["flavor_operation"](table[n])
        def entry(i, j):
            return reduce_rational(v[i,j])+256*reduce_rational(v[i+3,j])
        # Use the first two coordinates on the plane x+5y+29z=0.
        tangent = [(entry(i,j)-entry(i,2)*(1,5)[j]*pow(29,-1,prime)) % prime
                   for i in range(2) for j in range(2)]
        finite_pairs.append([x % prime for x in pair_coordinates(tangent)])
    determinant = abs(acb_mat(analytic_pairs).det())
    residue = int(nmod_mat(finite_pairs, prime).det())
    assert determinant > arb("1.16e-50") and residue != 0
    print("Collective pair determinant",determinant.str(10),"finite residue",residue)


if __name__ == "__main__":
    verify()

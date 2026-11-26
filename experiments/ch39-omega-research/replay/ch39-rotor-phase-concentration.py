"""Finite-rotor concentration at the Voronoi cell of the 10875th zeta phase."""
from flint import acb, arb, arb_mat, ctx


def window(bits):
    total = 2**bits
    angles = [2*(1/(2*z.imag)).atan() for z in acb.zeta_zeros(10874,3)]
    lower = (angles[1]+angles[2])*total/(4*arb.pi())
    upper = (angles[0]+angles[1])*total/(4*arb.pi())
    return list(range(int(lower.ceil().unique_fmpz()), int(upper.floor().unique_fmpz())+1))


def positive_definite(matrix):
    for size in range(1, matrix.nrows()+1):
        det = arb_mat([[matrix[i,j] for j in range(size)] for i in range(size)]).det()
        if det < 0:
            return False
        assert det > 0
    return True


def concentration_below(total, positions, cutoff, parity):
    def kernel(d):
        if d == 0:
            return arb(2*cutoff+1)/total
        theta = arb.pi()*d/total
        return ((2*cutoff+1)*theta).sin()/(total*theta.sin())
    return positive_definite(arb_mat([[
        (arb(".99") if i == j else 0)-kernel(i-j)-parity*kernel(i+j)
        for j in positions] for i in positions]))


def verify():
    ctx.prec = 160
    for bits, first, last, one, even, odd in [(32,64358,64363,765237620,765232862,765240547),
                                             (33,128716,128726,860203871,860199462,860204267)]:
        total = 2**bits
        positions = window(bits)
        assert positions == list(range(first,last+1))
        for parity,boundary in [(0,one),(1,even),(-1,odd)]:
            assert concentration_below(total,positions,boundary-1,parity)
            assert not concentration_below(total,positions,boundary,parity)
        exposure = min(even,odd)
        print("Paired-window threshold", bits, "bits, flux", exposure,
              "minimum closed-hop lower bound", 2*exposure+exposure%2)
    positive = window(32)
    cos_cap = max((2*arb.pi()*33366*j/(2**32)).cos().upper() for j in positive)
    assert 1/(1-cos_cap) < arb(".500000004")
    print("Nonnegative flux states: paired-window probability <=", (1/(1-cos_cap)).upper().str(15))


if __name__ == "__main__":
    verify()

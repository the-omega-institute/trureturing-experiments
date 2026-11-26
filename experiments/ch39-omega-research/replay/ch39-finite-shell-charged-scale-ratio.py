"""One-loop optimum for the parity-marginal heavy completion and the shell doublets."""
from fractions import Fraction as Q

from flint import arb, ctx


def main():
    ctx.prec = 192
    # Three Dirac color triplets with |Y|=2/3 and one charged Dirac singlet.
    assert 3*Q(4,3)*Q(3,5)*3*Q(2,3)**2==Q(16,5)
    assert 3*Q(4,3)*Q(1,2)==2 and Q(4,3)*Q(3,5)==Q(4,5)
    h = arb(125)**2
    inverse = [8*arb.pi()**2/arb(v)**2 for v in (".467453",".63811","1.0583")]
    k = 5*inverse[0]-3*inverse[1]-2*inverse[2]
    f = arb(109)/15*(inverse[1]-inverse[2])-arb(23)/6*(inverse[0]-inverse[1])
    assert f>254 and k<1084
    a,b = arb(8).log()/2,arb(24).log()/2
    least = 3-arb(5).sqrt()
    assert least>arb(3)/4 and -least.log()/2<arb(1)/4 and b<arb("1.6")
    # If the heaviest charged fermion has not opened, the 23-doublet
    # exposure alone forces a larger mass ratio than the interior optimum.
    assert (15*f/(22*23)-arb(1)/4-b).exp()>280
    ms = ((1000**2-h)/least).sqrt()
    assert ms>arb("1135.14915") and ms<arb("1135.14916")
    eigenvalues = [arb(6)]*7+[arb(4)]*6
    for j,n in ((5,3),(3,2)):
        for sign in (-1,1):
            eigenvalues += [3+sign*arb(j).sqrt()]*n
    z = sum((value+h/ms**2).log() for value in eigenvalues)/2
    A = k-44*(ms/1000).log()+12*a+4*b
    B = 15*f+22*z-402*a-46*b
    det = 44*448+16*506
    ell,w = (448*A+16*B)/det,(506*A-44*B)/det
    ratio = (ell-w).exp()
    crossing = ms*ell.exp()
    mf = ms*ratio
    assert ratio>arb("266.6444299") and ratio<arb("266.6444300")
    assert mf>arb("302681.1992") and mf<arb("302681.1994")
    assert crossing>440491221395 and crossing<440491221398
    assert w>b and ell>arb(6).log()/2+1
    assert Q(58*44)-60*22*Q(23,64)>0
    # Check all three inverse couplings independently at the attained endpoint.
    beta = (arb(41)/10,-arb(19)/6,-arb(7))
    scalar = (arb(1)/10,arb(1)/6,arb(0))
    color = (arb(16)/5,arb(0),arb(2))
    singlet = (arb(4)/5,arb(0),arb(0))
    total = sum((crossing/(h+val*ms**2).sqrt()).log() for val in eigenvalues)
    result = [inverse[i]-beta[i]*(crossing/1000).log()-scalar[i]*total
              -color[i]*(crossing/(mf*arb(8).sqrt())).log()
              -singlet[i]*(crossing/(mf*arb(24).sqrt())).log() for i in range(3)]
    assert all(abs(result[i]-result[j])<arb("1e-48") for i,j in ((0,1),(0,2),(1,2)))
    print("Minimum MF/MS",ratio.str(15),"at scalar scale",ms.str(14),
          "GeV, fermion scale",mf.str(14),"GeV, crossing",crossing.str(16),"GeV.")


if __name__=="__main__":
    main()

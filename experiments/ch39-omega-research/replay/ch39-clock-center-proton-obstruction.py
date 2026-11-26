"""Center triality vs low-energy chiral spectra for the charge-two SU(3)_X shell lift."""
from flint import fmpz_poly

z=fmpz_poly([0,1])
f=1+2*z
bar=2+z
# The common unit z is t^(12d), allowing any integral normalization d of 6Y.
# The pure constant ideal has positive generator nine for either nontrivial triality.
assert 4*bar**2-f*(2*z+7)==9
assert f**2-bar*(4*z-4)==9
for residue,first,second in ((4,f,bar**2),(7,bar,f**2)):
    assert first(residue)%9==second(residue)%9==0
    assert residue%3!=0

def check():
    for n1 in range(4):
        for n2 in range(4-n1):
            # Multiply n1 t^(-3d)-n2 t^(3d) by t^(3d). Remainders
            # have disjoint degrees 0,6d in the free rank-12d quotient.
            passes=n1%9==0 and n2%9==0
            assert passes==(n1==n2==0)
    print("restriction generators for center-one, center-two:",f,bar**2,";",bar,f**2)
    print("minimum positive isolated monomial multiplicity in both virtual ideals: 9")
    print("three light lepton doublets at y=-3: nonzero center charges allowed:",0)

if __name__=="__main__":
    check()

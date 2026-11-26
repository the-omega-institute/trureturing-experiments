"""A representation-ring parity obstruction for decoupling all bulk colored weak exotics."""
from dataclasses import dataclass
from itertools import product

@dataclass(frozen=True)
class Rep:
    dim: int
    boxes: int
    index: int
    cubic: int

    def tensor(self,o):
        return Rep(self.dim*o.dim,(self.boxes+o.boxes)%4,self.index*o.dim+self.dim*o.index,
                   self.cubic*o.dim+self.dim*o.cubic)

    def check(self):
        assert self.dim*self.boxes%4==0
        assert (self.index-self.boxes*self.cubic)%2==0


def main():
    # Traces on a traceless diagonal H, evaluated with all weights of the six.
    h=(1,1,1,-3)
    six=[h[i]+h[j] for i in range(4) for j in range(i+1,4)]
    assert sum(x**2 for x in six)==2*sum(x**2 for x in h)
    assert sum(x**3 for x in six)==0 and sum(x**3 for x in h)!=0
    generators=[Rep(4,1,1,1),Rep(6,2,2,0),Rep(4,3,1,-1)]
    for g in generators:
        g.check()
    # The ring proof extends to all tensor words: only the following two universal
    # cross-terms are left after using additivity of trace indices for tensors.
    for R,S in product(generators,repeat=2):
        assert (R.boxes*R.dim*S.cubic+S.boxes*S.dim*R.cubic)%2==0
        R.tensor(S).check()
    # Weak index 2T for the d-dimensional SU(2) irrep, normalized to one on a doublet.
    assert (2*(2*2-1))//6==1
    # Original and extra SU4 fundamentals of the explicit anomaly-free example.
    signatures=(3,9)
    assert sum(signatures)%2==0 and signatures[0]%2==1
    print("SU4 generators (boxes,cubic index,quadratic index)",[(r.boxes,r.cubic,r.index) for r in generators])
    print("three-family light color-weak signature",signatures[0],"extra signature in example",signatures[1])


if __name__=="__main__":
    main()

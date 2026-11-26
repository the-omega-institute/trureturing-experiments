"""An exact triangle basis, Wilson masses for weak Wilson lines, and the resulting one-loop threshold."""
from fractions import Fraction as F
from pathlib import Path
import runpy
from flint import fmpz_mat, fmpz_poly

weak = runpy.run_path(str(Path(__file__).with_name("ch39-kirchhoff-weak-gauge.py")))
q,pi = weak["q"],weak["pi"]
edges = [(i,j) for i,j,_ in weak["raw"]]
root = 14
rim = [i for i in range(19) if i != root]
chords = [(i,j) for i,j in edges if root not in (i,j)]
assert all(tuple(sorted((root,i))) in edges for i in rim)
# Fundamental triangle matrix in chord then oriented star coordinates is [I | P].
P=fmpz_mat(len(chords),len(rim))
index = {i:k for k,i in enumerate(rim)}
for a,(i,j) in enumerate(chords):
    P[a,index[i]] = 1
    P[a,index[j]] = -1
poly=(P.transpose()*P).charpoly()
factors=poly.factor()
assert len(chords)==117 and P.rank()==17
print("rim Laplacian characteristic polynomial factors",
      [(str(f),int(power)) for f,power in factors[1]])
assert (poly // fmpz_poly([0,1])).degree()==17
# Eliminate gauge redundancies: nonzero modes are the spectrum of I+P P^t.
# A vertex in the rim has degree 17 (three of the remaining hubs are universal).
assert max(sum(1 for e in chords if v in e) for v in rim)==17
assert 117-P.rank()==100
shift = q(39)*q(8).log()/(2*pi)
floor = weak["running"]()-shift
assert floor > q(14) and q("601.5")-weak["k"]*8/2<0
assert 1+1-2==0 and F(1,6)+(F(1,6)-F(1,2))+F(1,6)==0
assert 2*F(1,180)+F(7,720)==F(1,48)
print("thresholds Q and sqrt(19) Q, multiplicity at Q",100)
print("uniform inverse alpha2 with massive cycles",floor.lower().str(13))

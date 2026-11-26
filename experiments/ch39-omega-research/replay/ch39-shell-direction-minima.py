"""Which positive direction vacua survive the clock ordering, and how do their rotational curvatures differ?"""
from pathlib import Path
from flint import arb,fmpz_mat
import runpy

m=runpy.run_path(str(Path(__file__).with_name("ch39-shell-vacuum-directions.py")))
X,T=m["X"],m["T"]
g=m["g"]
def sq(x,y):
    return sum((a-b)**2 for a,b in zip(x,y))
reps={tuple((tuple(sorted(seq[j])) for j in (0,6))):seq for seq in m["solutions"] if
      (seq[0][0]<seq[0][1]<seq[0][2])}
assert len(reps)==4
found=[]
graphs=[]
for shape,Y in sorted(reps.items()):
    C=fmpz_mat([[sum(y[i]*x[j] for x,y in zip(X,Y)) for j in range(3)] for i in range(3)])
    a,b=C[0,0]-C[0,1],C[0,0]+2*C[0,1]
    den=g["den"]
    # This sign convention is fixed by the order of the three clock chords on the source shell.
    chords=[sq(Y[0],Y[j]) for j in (6,7,8)]
    assert sorted(sq(X[0],X[j]) for j in (6,7,8))==[6,12,14]
    source_order=[(sq(X[0],X[j]),sq(Y[0],Y[j])) for j in (6,7,8)]
    source_order.sort()
    candidates=[v for k,v in source_order]
    # Test whether the maximum arc allowed by the 19-shell graph preserves the six outer-to-inner
    # incidences for each outer basepoint as in the chosen common subgraph.
    image_chords={sq(Y[i],Y[j]) for i in range(6) for j in range(6,9)}
    mismatches=[(i,j) for i in range(9) for j in range(i)
                if (sq(X[i],X[j])<=14)!=(sq(Y[i],Y[j])<=804)]
    print("induced adjacency disagreements",shape,len(mismatches))
    if not mismatches:
        graphs.append(shape)
    if max(image_chords)<=804:
        found.append(shape)
    print(shape,"transverse / longitudinal increments",a,b,"three connecting chords",candidates,
          "oscillator sum without common scale",(2*(arb(int(a+b))/den).sqrt()+
          (2*arb(int(a))/den).sqrt()).str(11))
print("zero-direction vacua with all outer-inner clock arcs retained",found)

assert graphs==[((1,5,29),(13,13,23))]

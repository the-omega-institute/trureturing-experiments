"""Weak scalar spectral weights and the open decay into the two lowest vector harmonics."""
from pathlib import Path
import runpy
from flint import arb

m=runpy.run_path(str(Path(__file__).with_name("ch39-clock-weak-cubic.py")))
q,pi=m["q"],arb.pi()
z1,z2,z3=m["roots"]()
s=m["s"]
amp=q("-1.09974").union(q("-1.01335"))
p=q(".00617").union(q(".19955"))
WJ=q(".3337").union(q("6.430"))
wv=q("1.19576").union(q("1.21744"))

def weight(z):
    r=s*z/((z-1)*(z-11)-s*(2*z-11))
    # Partial fraction residues of the unperturbed coupling vectors' Gram resolvent.
    norm=(q(".9")*(r*r+1)-q(".2")*r)+q("1.1")*(r+1)**2*((z-1)/(z-11))**2
    v=r+amp
    return (p.lower()*v.abs_lower()**2/norm.upper()).lower().union(
           (p.upper()*v.abs_upper()**2/norm.lower()).upper())

def upper_inverse(x):
    # Positive harmonics of the no-vertex-kinetic-term weak graph interlace
    # its rank-one Dirichlet levels with at most 18 additional positive modes.
    shifts=q(0)
    for _,_,L in m["cur"]["weak"]["raw"]:
        first=2*pi/L
        for n in range(1,4):
            if not first*n>x:
                shifts+=max(q(0),(x/(first*n)).log().upper())
        assert 4*first>x
    # Eighteen right modes are already active; the extra cycle scalars screen.
    return q("40.9906")+(q("63.5")*x.log()+4*shifts)/(2*pi)

def check():
    lines=[weight(z) for z in (z1,z2,z3)]
    assert all(f>0 for f in lines)
    assert z2.sqrt()<2*wv<z3.sqrt()
    mass=z3.sqrt()
    inv=upper_inverse(mass)
    gs=(4*pi/inv.upper()).lower().union((4*pi/q("14.0866")).upper())
    gmin=gs.lower()*WJ.lower()*lines[2].lower()
    gmax=gs.upper()*WJ.upper()*lines[2].upper()
    spin=2+((z3-2*wv**2)/(2*wv**2))**2
    beta=(1-4*wv**2/z3).sqrt()
    # Both distinct adjoint color assignments for one fixed scalar color.
    k=spin*beta/(8*pi*mass)
    width=(gmin*k.lower()).lower().union((gmax*k.upper()).upper())
    print("bright weights",[(float(x.lower()),float(x.upper())) for x in lines])
    print("inverse alpha_2 upper at heavy line",inv.upper().str(12))
    print("partial width / Q",width.lower().str(12),width.upper().str(12))
    print("M Gamma / Q^2",(mass.lower()*width.lower()).lower().str(12),(mass.upper()*width.upper()).upper().str(12))

if __name__=="__main__":
    check()

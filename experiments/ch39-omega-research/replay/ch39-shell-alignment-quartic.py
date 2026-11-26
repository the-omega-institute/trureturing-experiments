"""Exact diagonal self-coupling and unsubtracted vacuum energy of the alignment action."""
from pathlib import Path
import runpy

newton=runpy.run_path(str(Path(__file__).with_name("ch39-shell-alignment-newton.py")))
q,pi,nu,Q=newton["q"],newton["pi"],newton["nu"],newton["Q"]
den=q(17*867).sqrt()
E0=18-2178/den
relative=nu*den/198
# A unit-generator motion with canonical coordinate phi=f theta has
# V=K E0+m^2 f^2*(1-cos(phi/f)), hence V''''(0)=-m^2/f^2.
fmin=3*nu.sqrt()/(4*pi).sqrt()
Kmin=(relative.lower()*fmin.lower()**2).union(relative.upper()*fmin.upper()**2)
absolute=(E0*Kmin.lower()*Q.lower()**4).union(E0*Kmin.upper()*Q.upper()**4)
for label,z in (("f/Q at quartic=4pi for the heavy pair",fmin),
                ("minimum K/Q^4 at that bound",Kmin),
                ("minimum unsubtracted alignment energy in Q^4",E0*Kmin),
                ("minimum unsubtracted alignment energy in GeV^4",absolute)):
    print(label,z.lower().str(13),z.upper().str(13))
assert fmin>q(199) and E0*Kmin>q("8e7")

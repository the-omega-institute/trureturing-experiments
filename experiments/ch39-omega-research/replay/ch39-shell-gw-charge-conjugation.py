"""Local charge-conjugation projection and the interval of Wilson masses at both extrema."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path
from flint import acb, acb_mat, arb, ctx


s=run_path(str(Path(__file__).with_name("ch39-shell-gw-symmetric-links.py")))
ctx.prec=240
model=s["model"]
data=model["data"]
pts,faces=data["maps"]()
sites,H,_=data["shell"]["make_operator"]()
one=arb(1); h=(one/2).sqrt()
spin=((one,0),(0,one),(h,h),(h,-h),(h,acb(0,h)),(h,acb(0,-h)))
neutral=[]
for i in range(12):
    row=[]
    for j in range(12):
        a,b=map(acb,spin[j%6]); c,d=map(acb,spin[i%6])
        phase=a.conjugate()*c+b.conjugate()*d
        value=data["ball"](H[i,j])
        row.append(value*phase.conjugate()/abs(phase) if H[i,j]!=0 else acb(0))
    neutral.append(row)
charge_symmetric=[[(neutral[i][j]-neutral[(i+6)%12][(j+6)%12].conjugate())/2 for j in range(12)]
                  for i in range(12)]


def scalar_kernel(sign,low,berry=True):
    flux=[]
    for face in faces:
        ds=[26-sum(a*b for a,b in zip(pts[i],pts[j])) for i,j in data["edges"](face)]
        flux.append(sign*(F(1,4) if ds==[7]*3 else F(5,12)))
    field=data["cochain"](faces,flux)
    k=[]
    for i in range(12):
        row=[]
        for j in range(12):
            p,v=pts[j%6],pts[i%6]
            x,y,z=p;a,b,c=v
            turn=field.get((j%6,i%6),F(0))
            phase=acb(0,2*arb.pi()*turn.numerator/turn.denominator).exp()
            if berry:
                ov=acb(26+z*c+x*a+y*b+arb(26).sqrt()*(z+c),x*b-y*a)
                phase*=ov/abs(ov)
            if low and 26-x*a-y*b-z*c==1 and p.index(1)!=0:
                phase=-phase
            row.append(phase*charge_symmetric[i][j])
        k.append(row)
    return k


if __name__=="__main__":
    for low in (False,True):
        for sign in (1,-1,0):
            k=scalar_kernel(sign,low)
            poly=acb_mat([[(-1 if i<6 else 1)*v for v in row] for i,row in enumerate(k)]).charpoly()
            roots=poly.roots(tol=1e-14)
            real=sorted([z.real for z in roots if z.imag.contains(0)],key=lambda x:float(x.mid()))
            assert len(real)==(2 if sign else 0)
            if sign:
                assert 0<real[0]<arb(".16") and arb(".63")<real[1]<2
            for i in range(12):
                k[i][i]+=one/2 if i<6 else -one/2
            gap=one/200
            assert model["count"](k,-gap)==model["count"](k,gap)==6-sign
            s["calc"]["opposite_modes"](k,sign,gap)
            print("Wilson minimum" if low else "untwisted",sign,"index",sign,
                  "mass walls",[r.str(13) for r in real])

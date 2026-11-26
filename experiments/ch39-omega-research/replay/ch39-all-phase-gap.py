"""Uniform gap for arbitrary pairs of unitary two-component clocks, via three complex vectors."""
from pathlib import Path
from flint import arb,ctx,acb,acb_mat
root=Path(__file__).parent
code=(root/'ch39-quaternion-spectral-polynomial.py').read_text()
code=code[:code.index('def evaluate(')]
code=code.replace('out=one.copy()','out=[]')
a=code.index('        re=times(one,c0.real)')
code=code[:a]+'''        out.append((c0.real,c))
    return H, block(O,O), out
'''
source={'__file__':str(root/'ch39-quaternion-spectral-polynomial.py'),'__name__':'sector'}
exec(compile(code,'sectors','exec'),source)
src=source['src']
def gap_at(d,energy,steps=256):
    src['rho']=arb.pi()**2/(d*src['theta'](804))
    H,A,blocks=source['polynomial'](-energy)
    # Positivity of the outside principal block is a separate condition.
    for j in range(1,A.nrows()+1):
        det=acb_mat([[A[a,b] for b in range(j)] for a in range(j)]).det()
        if not det.real>0:return None
    estimates=[]
    for const,(c1,c2,c3) in blocks:
        low=[]
        for j in range(steps):
            p=(2*j+1)*arb.pi()/steps
            z=acb(0,p).exp()
            value=const+2*((c1*z).real-abs(c2+c3*z.conjugate()))
            # Use the Lipschitz derivative of both terms in the remaining angle.
            low.append((value-2*(abs(c1)+abs(c3))*arb.pi()/steps).lower())
        estimates.append(min(low))
    return min(estimates)
if __name__=='__main__':
    for energy in ('.3','.45','.6','.75','.90','1.0','1.1'):
        ans=gap_at(arb(256),arb(energy))
        print(energy,ans.str(10) if ans is not None else 'outside indefinite')

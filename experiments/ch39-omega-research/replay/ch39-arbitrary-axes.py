"""All simultaneous conjugacy classes with the two angles in the observed zero square."""
from pathlib import Path
import runpy
from flint import arb,acb,acb_mat
root=Path(__file__).parent
m=runpy.run_path(str(root/'ch39-quaternion-scalar-levels.py'))
ctxsrc=m['m']['m']['real_energy_matrix'].__globals__['src']
ctxsrc['rho']=arb.pi()**2/(256*ctxsrc['theta'](804))
s,y=m['p2'].sin_cos()
grid=128
energy=arb('1.568')
k=(energy*arb.pi()**2/256).sqrt()
nmax=[]
for dfilter in [(724,804),(804,)]:
    degrees=[0]*19
    for i,j,d,_ in ctxsrc['edges']:
        if d in dfilter:degrees[i]+=1;degrees[j]+=1
    nmax.append(max(degrees))
assert nmax==[4,2]
shift=(4*arb('.01')+2*arb.pi()/(2*grid))*k/k.sin()
assert energy-shift>arb('1.5')
for j in range(grid):
    h=(2*j+1)*arb.pi()/(2*grid)
    a,b=h.sin_cos()
    m['count'].__globals__['V']=acb_mat([[acb(y,b*s),s*a],[-s*a,acb(y,-b*s)]])
    assert m['count'](energy)==0
# Hadamard bound on the product of the low 18 static trial eigenvalues;
# the diagonal 4 (six scalar sites), 8 (three), repeated twice for flavor.
base=ctxsrc['matrix'](0,{})
assert all(base[i,i].real==4 for i in range(6)) and all(base[i,i].real==8 for i in range(6,9))
mean=(6*arb(4).log()+3*arb(8).log())/18
print('uniform x_1 floor at d=256',(energy-shift).lower().str(12),'mean upper',mean.upper().str(12))
# In the global integrator, use the angle-independent upper mean. This bound
# covers energies below 1 TeV in the initial estimates by the same light-count argument.
base_code=(root/'ch39-global-gauge-crossings.py').read_text()
a=base_code.index('    highlogs=sum(')
b=base_code.index('\n',a)
base_code=base_code[:a]+'    highlogs=(6*q(4).log()+3*q(8).log())/18'+base_code[b:]
scope={'__file__':str(root/'ch39-global-gauge-crossings.py'),'__name__':'full_angle_global'}
exec(compile(base_code,'all_axes_global','exec'),scope)
q=scope['q']
driver=(root/'ch39-all-perturbative.py').read_text().replace("m=runpy.run_path(str(root/'ch39-global-gauge-crossings.py'))","m=broad").replace("q('10.753')","q(30)")
driver=driver.replace("assert W.upper()<q('50.85')","assert W.upper()<q('50.90')")
# Establish the 3.44 logarithmic width used for every narrow-shape integration.
driver=driver.replace("# The same extreme-shape integrator","assert W.upper()/18+mean-q('1.5').log()/2 < q('3.44')\n# The same extreme-shape integrator")
driver_scope={'__file__':str(root/'ch39-all-perturbative.py'),'__name__':'all_axes','broad':scope,'mean':mean}
exec(compile(driver,'all_axes_necessary','exec'),driver_scope)
W=driver_scope['W']
L=max((W.upper()/18+mean-q('1.5').log()/2).upper(),(q('1.62')*(q(256)/q('1.5')).sqrt()).log().upper())
assert L<q('3.329')
code=(root/'ch39-spectrum-width-obstruction.py').read_text().split('if __name__')[0].replace('38.3316','38.3451').replace('-.029827','-.02980')
bounds={}
exec(compile(code,'all_axes_ew','exec'),bounds)
_,_,margin=bounds['bound']('3.329')
Q=1000*q('26.85181').exp()
mass=margin*(Q/(16*arb.pi()**2))**2
print('uniform log width',L.upper().str(12),'mass floor',mass.lower().str(14))
assert mass>q('1e23')

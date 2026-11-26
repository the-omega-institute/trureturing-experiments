"""All clock phases and all relative axes: conditional perturbative global bounds plus crossing existence."""
from pathlib import Path
import runpy
from flint import arb
root=Path(__file__).parent
sector=runpy.run_path(str(root/'ch39-all-phase-gap.py'))
q=lambda s: arb(s)
dstar=q(384);xmin=q('.92')
gap=sector['gap_at'](dstar,xmin,1024)
assert gap>0
print('three-correlation gap',gap.str(12))
mean=(6*q(4).log()+3*q(8).log())/18
base_code=(root/'ch39-global-gauge-crossings.py').read_text()
a=base_code.index('    highlogs=sum(');b=base_code.index('\n',a)
base_code=base_code[:a]+'    highlogs=(6*q(4).log()+3*q(8).log())/18'+base_code[b:]
base_code=base_code.replace("q('1.5')","q('.92')").replace('q(256)','q(384)')
scope={'__file__':str(root/'ch39-global-gauge-crossings.py'),'__name__':'all_phases_coarse'}
exec(compile(base_code,'all_phases_coarse','exec'),scope)
driver=(root/'ch39-all-perturbative.py').read_text().replace("m=runpy.run_path(str(root/'ch39-global-gauge-crossings.py'))","m=broad")
driver=driver.replace("q('10.753')","q(30)").replace("q('1.5')","q('.92')").replace('q(256)','q(384)').replace('14**2','22**2').replace('q(14).log()','q(22).log()')
driver=driver.replace("assert W.upper()<q('50.85')","assert W.upper()<q('51.1')")
driver=driver.replace("env=env[:env.index('def certify(')]","env=env[:env.index('def certify(')].replace(\"width=q('3.44')\",\"width=q('3.71')\")")
driver=driver.replace("# The same extreme-shape integrator","assert W.upper()/18+mean-q('.92').log()/2<q('3.71')\nassert (q('1.62')*(q(384)/q('.92')).sqrt()).log()<q('3.71')\n# The same extreme-shape integrator")
check="assert T>q('26.837') and T<q('26.899') and W<q('49.222') and alpha<q('38.364')"
driver=driver.replace(check,"assert T>q('26.845') and T<q('26.889') and W<q('49.2') and alpha<q('38.36')")
h={'__file__':str(root/'ch39-all-perturbative.py'),'__name__':'all_phases_main','broad':scope,'mean':mean}
exec(compile(driver,'all_phases_main','exec'),h)
W=h['W']
assert W>q('48.23') and W<q('49.02') and h['alpha']>q(38)
r=(W.upper()/18+mean).exp()
L=max((r/xmin.sqrt()).log().upper(),(q('1.62')*(dstar/xmin).sqrt()).log().upper())
assert r<35 and L<q('3.59')
# For two fermion modes, t=min exposures must be >=0.02. Otherwise
# the second level of the lightest line contributes more than 0.02.
assert 35*(-q('.98')).exp()<q('2.22')*(1+q('2.22')**2)
assert q(2).sqrt()*(1+2*q('2.22')**2)/(1+q('2.22')**2)<q('.96').exp()
def wt(t):return (1+12*t)*(-2*t).exp()
assert wt(q('.02'))+wt(q('.98'))<2+13*(-q(2)).exp()
code=(root/'ch39-spectrum-width-obstruction.py').read_text().split("if __name__")[0]
code=code.replace('38.3316','38.36').replace('-.029827','-.0295')
code=code.replace('w2=1+13*(-q(2)).exp()', "w2=(1+12*q('.02'))*(-q('.04')).exp()+(1+12*q('.98'))*(-q('1.96')).exp()")
bounds={}
exec(compile(code,'all_phases_width','exec'),bounds)
coef=bounds['bound']('3.59')[2]
Q=1000*q('26.845').exp()
mass=coef*(Q/(16*arb.pi()**2))**2
assert mass>q('1e23')
print('r, width, dimensionless and physical floors',r.upper().str(12),L.upper().str(12),coef.lower().str(12),mass.lower().str(14))
# All continuous threshold shapes satisfy the same broad rectangle existence test.
env=(root/'ch39-exposure-envelope.py').read_text().split("if __name__")[0].replace("width=q('3.44')","width=q('3.71')")
e={'__file__':str(root/'ch39-exposure-envelope.py')}
exec(compile(env,'all_phases_existence','exec'),e)
e['certify']('21.425','21.52','26.837','26.899')

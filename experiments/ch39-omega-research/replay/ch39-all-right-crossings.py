"""Perturbative r=400 crossings, allowing arbitrary weak-singlet right threshold profiles."""
from fractions import Fraction
from pathlib import Path
import runpy

from flint import arb

root = Path(__file__).resolve().parent
model = runpy.run_path(str(root/'ch39-right-scalar-rge.py'))
m = model['data']['g']
g = model['g']
q, D, cap = g['q'], g['loop'], m['cap']
L = model['Ls']
dot, cut = m['dot'], m['cut']
two, four = q(2), q(4)
def coefficient(i, counts):
    nl, nr, nf = counts
    return [a+nl*l+nr*r+nf*f for a,l,r,f in zip(g['B'][i],g['BD'][i],model['BR'][i],g['BF'][i])]
def value(row, gs):
    return dot(row,gs)/D
# First-exit weak bound. Extra right states have no weak row.
trial = [cap,two,cap]
assert g['beta'][1]+value(g['B'][1],trial) < 0
loss = (g['d'][1]+value(g['BD'][1],trial))*model['WL']+(8+value(g['BF'][1],trial))/3
assert g['G0'][1]-2*loss/D > 1/two
assert g['beta'][0]-sum(g['D'][0])*3*cap/D > 0
# At first hitting 1.4 the color right derivative must be nonnegative.
# If it later hits 4, its derivative stays positive there for all subsequent counts.
def color_rates(nl,nr,nf, first, second, ew, yuk):
    row = coefficient(2,(nl,nr,nf))
    assert row[0] >= 0 and row[1] >= 0
    leading = -7+(q(nl)+nr)/3+8*nf
    upper = leading+value(row,[ew,ew,first])
    lower = leading+(row[2]*second-sum(g['D'][2])*yuk)/D
    return upper, lower
def bound_color(first, second, ew, yuk):
    assert 1/g['G0'][2] < first and ew < second
    for nl in range(19):
        for nf in (q(0),q('1/3')):
            for nr in range(32):
                upper, lower = color_rates(nl,nr,nf,first,second,ew,yuk)
                assert upper < 0 or lower > 0
                if nr==31: assert lower > 0
bound_color(q('1.4'),four,two,3*cap)

# At a crossing the projected integrated equations give tau and WL-WR.
# Move the averaged right term in the latter equation to the denominator.
rcR = [m['proj'](row,model['BR']) for row in m['P']]
beta_floor = g['beta'][0]-sum(g['D'][0])*3*cap/D
initial_time = m['nonneg'](g['G0'][0]*D/(2*beta_floor))
initial_right = m['nonneg'](g['G0'][0]*D/(2*q('8/15')))
initial_gauge = [m['nonneg'](two),m['nonneg'](two),m['nonneg'](four)]
initial_yukawa = [m['nonneg'](3*cap)]*3

def intersect_flow(tau,right,gs,ytr,steps,left_logs=L):
    def smallest(a,b):
        return min(a.lower(),b.lower()).union(min(a.upper(),b.upper()))
    wl = sum(smallest(tau,li) for li in left_logs)
    wf = smallest(tau,q(1))/3
    gs_avg = [q(0)]*3
    union_g = union_l = union_f = None
    for j in range(steps):
        frac = q(Fraction(j,steps)).union(q(Fraction(j+1,steps)))
        dtau = tau*frac
        rel = tau*(frac-1)
        el = m['nonneg'](min((frac*wl).upper(),sum(g['positive'](li+rel) for li in left_logs).upper()))
        er = m['nonneg']((frac*right).upper())
        ef = m['nonneg'](g['positive'](1+rel).upper()/3)
        if tau.lower() > max(li.upper() for li in left_logs):
            el = sum(g['positive'](li+rel) for li in left_logs)
            ef = g['positive'](1+rel)/3
        box = []
        for i,ri in enumerate(map(q,('8/15',0,'1/3'))):
            sm = dot(g['B'][i],gs)-dot(g['D'][i],ytr)
            heavy = dot(g['BD'][i],gs)*el+dot(model['BR'][i],gs)*er+dot(g['BF'][i],gs)*ef
            one = g['beta'][i]*dtau+g['d'][i]*el+ri*er+8*ef
            G = g['G0'][i]-2*(one+(sm*dtau+heavy)/D)/D
            low = max(1/gs[i].upper(),G.lower())
            assert G.upper() >= low
            box.append(cut(gs[i],1/low.union(G.upper())))
        gs_avg = [a+tau*b/steps for a,b in zip(gs_avg,box)]
        def merge(old):
            return box if old is None else [x.union(y) for x,y in zip(old,box)]
        union_g = merge(union_g)
        if not rel+max(li.upper() for li in left_logs)<0: union_l = merge(union_l)
        if not rel+1<0: union_f = merge(union_f)
    nonright = m['induced'](tau,wl,gs_avg,[tau*y for y in ytr],[wl*z for z in union_l],[wf*z for z in union_f])
    R = cut(right,(wl-nonright[1])/(1-dot(rcR[1],union_g)/D))
    T = cut(tau,nonright[0]-R*dot(rcR[0],union_g)/D)
    return T,R,union_g

def crossings(left_logs):
    T,W,gs,ys = initial_time,initial_right,initial_gauge,initial_yukawa
    for steps in (8,32,128,256,256):
        T,W,gs = intersect_flow(T,W,gs,ys,steps,left_logs)
    # One-loop logarithmic light Yukawas, using all three earlier gauge floors.
    assert all(lo>q(0) for lo in (v.lower() for v in gs))
    y = [z.upper() for z in g['y0']]
    ytop = [q(0)]*3
    step = T.upper()/256
    for j in range(256):
        tube = [(v+q('.01')).exp() for v in y]
        tri = [sum(tube[i:i+3]) for i in (0,3,6)]
        total = 3*(tri[0]+tri[1])+tri[2]
        t = (j*step).union((j+1)*step)
        gauge_floor = []
        for i in range(3):
            adverse = sum(min(q(0),v) for v in g['B'][i])*cap-3*cap*sum(g['D'][i])
            upper_G = g['G0'][i]-2*t*(g['beta'][i]+adverse/D)/D
            gauge_floor.append(max(gs[i].lower(),1/upper_G.upper()))
        for i in range(9):
            rate = (3*tube[i]+2*(total-dot(g['C'][i//3],gauge_floor)))/D
            assert step*rate<q('.01')
            y[i]=(y[i]+step*rate).upper()
        ytop=[max(a,b.upper()) for a,b in zip(ytop,tri)]
    ys=[m['nonneg'](v) for v in ytop]
    for steps in (256,512,512):
        T,W,gs = intersect_flow(T,W,gs,ys,steps,left_logs)
    return T,W

T,W = crossings(L)
print('all crossing time and exposure',T.str(12),W.str(12))
assert T > q('26.89357') and T < q('27.08326')
assert W > q('42.47815') and W < q('47.80191')
read = runpy.run_path(str(root/'ch39-arbitrary-right-crossing.py'))
price, needed = read['required_exposure'](16,T.lower())
gap = (needed-W)*price
print('all perturbative crossings, u=16 necessary right logs above',needed.lower().str(13),
      'trace gap at Q/4 above',gap.lower().str(13))
assert needed > W and gap > 0

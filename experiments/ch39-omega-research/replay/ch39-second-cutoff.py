"""A second cutoff for the full or any truncated charged left tower at the r=400 crossing."""
from pathlib import Path
import runpy
from flint import arb

root = Path(__file__).resolve().parent
model = runpy.run_path(str(root/'ch39-finite-charge.py'))
q = model['q']
m = model['model']
factors = m['factors']
low = model['low_y']
h = m['h']
pi = arb.pi()
scale = q(400)**2
# Boundary-zero codimension 38, and 18 low states already accounted for.
def upper_high_heat():
    total = 20*(-q('1.2710')).exp()
    for mult, factor in factors:
        c = factor/scale
        for j in range(1,6):
            eigen_floor = c*j*j/(1+pi**2*j*j*h*h/6)
            total += mult*(-eigen_floor).exp()
        # For the remaining harmonics use y >= c*j^2/3 and j^2 >= 6*j.
        a = (-2*c).exp()
        total += mult*a**6/(1-a)
    return total

u = q(6)
right_y = (-2*q('2.55').union(q('2.605'))).exp()
low_left = sum((-y).exp() for y in low)
delta = low_left - 18*(-right_y).exp()
high = upper_high_heat()
def low_gap(u):
    return sum((-u*y).exp() for y in low)-18*(-u*right_y).exp()
universal = low_gap(u) - (-(u-1)).exp()*(delta+high)
# For a two-mode selection w_j > 0 of total W, sum w_j^2 <= W^2.
band = model['selected_band']
w = delta+2*(-band).exp()
second = low_gap(q(2))+2*(-2*band).exp()-w**2
print('full high heat upper',high.upper().str(12))
print('S_L(6)-S_R(6) lower',universal.lower().str(12))
print('two selected modes S_L(2)-S_R(2) lower',second.lower().str(12))
assert universal > q('.1') and second > q('.09')

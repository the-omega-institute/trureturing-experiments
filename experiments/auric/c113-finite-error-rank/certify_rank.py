#!/usr/bin/env python3
"""Exact finite coefficient certificates for the unchanged C113 5 by 5 table.
Uses rational arithmetic only for inequalities; no numerical sign samples.
Outputs definitions, exact coefficients, and finite assertions as reusable result data.
"""
import sys
assert sys.version_info >= (3, 12), sys.version
from fractions import Fraction as Q
from itertools import permutations, combinations
from math import factorial, prod
from pathlib import Path
import argparse, json, hashlib, subprocess
import sympy as sp

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--repository', type=Path, required=True,
                    help='Git repository containing the fixed C113 source commit')
parser.add_argument('--output', type=Path, required=True,
                    help='Destination JSON file; may be outside either repository')
args = parser.parse_args()
assert sp.__version__ == '1.14.0', sp.__version__
PIN = '2ecdfbb41b5bd234229ea9b61176a174a4c4de90'
SOURCE = 'docs/develop/theory/AURIC_INTERNAL_BALANCED_CALIBRATION_AND_FILTER_COMPOSITION_GEOMETRY_C113.md'
SOURCE_SHA256 = '39921aa9b28ea442cf0e0a20f93a6ff872207baec6738d1cb1cc5f3318e9db5a'
source_bytes = subprocess.check_output(['git', '-C', str(args.repository.resolve()),
                                      'show', f'{PIN}:{SOURCE}'])
assert len(source_bytes) == 163924
assert hashlib.sha256(source_bytes).hexdigest() == SOURCE_SHA256
N=7  # coefficients through t^7, exact quotient Q[t]/(t^8)
ZERO=(Q(0),)*(N+1)
ONE=(Q(1),)+(Q(0),)*N

def add(a,b): return tuple(x+y for x,y in zip(a,b))
def scale(a,c): return tuple(c*x for x in a)
def mul(a,b):
    r=[Q(0)]*(N+1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b[:N+1-i]):
                if y: r[i+j]+=x*y
    return tuple(r)
def det(a):
    n=len(a); ans=ZERO
    for p in permutations(range(n)):
        term=ONE
        for i in range(n): term=mul(term,a[i][p[i]])
        parity=sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))%2
        ans=add(ans,scale(term,(-1)**parity))
    return ans

def s(n):
    return tuple((Q(1,4) if k==0 else Q(0)) +
                 Q((-1)**k * (2**(2*k)+4**(2*k)) * (n+1)**(2*k),8*3**n*factorial(2*k))
                 for k in range(N+1))
H=[[s(i+j) for j in range(5)] for i in range(5)]
e={0:ONE}
for k in range(1,6):
    val=ZERO
    for I in combinations(range(5),k): val=add(val,det([[H[i][j] for j in I] for i in I]))
    e[k]=val
char={5-k:scale(e[k],(-1)**k) for k in range(6)}
F=det([row[:4] for row in H[:4]])
P=Q(2783,13122); A=Q(36515,1062882); B=-Q(50420,43046721); D=Q(256,43046721)
assert char[3][0]==P
assert char[2][0]==0 and char[2][1]==A
assert char[1][:3]==(0,0,0) and char[1][3]==B
assert char[0][:6]==(0,)*6 and char[0][6]==-D
assert F[:3]==(0,0,0) and F[3]==-Q(10,59049) and F[4]==-Q(136,177147)
assert -A/P==-Q(36515,225423)
assert -B/A==Q(20168,591543)
assert D/B==-Q(64,12605)

# Independent Newton identities check of every characteristic coefficient modulo t^8.
def matmul(a,b):
    return [[sum_polys([mul(a[i][k],b[k][j]) for k in range(5)]) for j in range(5)] for i in range(5)]
def sum_polys(xs):
    r=ZERO
    for x in xs:r=add(r,x)
    return r
power=H
tr={}
for k in range(1,6):
    tr[k]=sum_polys([power[i][i] for i in range(5)])
    if k<5:power=matmul(power,H)
for k in range(1,6):
    rhs=sum_polys([scale(mul(e[k-j],tr[j]),(-1)**(j-1)) for j in range(1,k+1)])
    assert scale(e[k],k)==rhs

# Exact generalized Vandermonde identities for every row subset of size 4.
y=sp.symbols('y0:4')
V=prod(y[j]-y[i] for i in range(4) for j in range(i+1,4))
schur={}
for missing in range(5):
    exps=[i for i in range(5) if i!=missing]
    m=4-missing
    elementary=sum(prod(y[i] for i in I) for I in combinations(range(4),m))
    M=sp.Matrix([[yy**a for yy in y] for a in exps])
    assert sp.Poly(sp.expand(M.det(method='domain-ge')-V*elementary),*y).is_zero
    schur[str(missing)]={'exponents':exps,'factor':f'e_{m}','bound_on_unit_polydisc':__import__('math').comb(4,m)}

# Exact modulus Vandermonde product in w=e^(2 i theta), conjugate w -> 1/w.
w=sp.symbols('w',nonzero=True)
z=[sp.Integer(1),w**(-2)/3,w**(-1)/3,w/3,w**2/3]
c=[sp.Rational(1,4),w**(-2)/16,w**(-1)/16,w/16,w**2/16]
V5=prod(z[j]-z[i] for i in range(5) for j in range(i+1,5))
assert sp.simplify(prod(c)-sp.Rational(1,2**18))==0
# |z_j-z_i|^2 is represented by d(w)*d(1/w) on |w|=1.
det_abs=prod(c)*V5*V5.subs(w,1/w)
sin2=lambda k:(2-w**k-w**(-k))/4  # sin^2(k theta)
q=lambda k:(10-3*(w**(k//2)+w**(-k//2)))/9
expected=sp.Rational(1,2**18)*sp.Rational(2,3)**12*sin2(2)*sin2(4)*sin2(1)**2*sin2(3)**2*q(2)**2*q(4)**2
assert sp.factor(det_abs-expected)==0
for n in range(9):
    moment=sp.Rational(1,4)+sp.Rational(1,16*3**n)*(w**(n+1)+w**(-n-1)+w**(2*n+2)+w**(-2*n-2))
    assert sp.expand(sum(cj*zj**n for cj,zj in zip(c,z))-moment)==0

# Every 4x4 H minor: enough terms to reproduce leading compound and F coefficients.
allminors={}
for r in range(5):
    for cidx in range(5):
        I=[i for i in range(5) if i!=r]; J=[j for j in range(5) if j!=cidx]
        d=det([[H[i][j] for j in J] for i in I])
        assert d[:3]==(0,0,0)
        allminors[f'delete_{r}_{cidx}']=[str(x) for x in d]

# Analytic remainder envelope inequalities at Tmax = 2^-26.
T=Q(1,2**26)
r4=Q(136,177147)*T/Q(10,59049)
r5=48*16**5*T**2/Q(10,59049)
assert 16*T<Q(1,2)
assert r4<Q(1,8) and r5<Q(1,8)
cluster_angles=[-2,-1,1,2]
triple_products=[prod(abs(a-b) for a,b in combinations(S,2)) for S in combinations(cluster_angles,3)]
assert max(triple_products)==12
assert prod(abs(a-b) for a,b in combinations([-4,-2,2,4],2))==4608
K=6*(Q(2048,243)+Q(512,81))
assert K==Q(7168,81)
minor_bound=5*K*K/Q(16384)
compound_bound=5*minor_bound
assert minor_bound==Q(15680,6561) and compound_bound==Q(78400,6561)
# det lower: sin x >= x/2 and q>=4/9; upper: sin x<=x, q<=16/9.
angle_leading=2**2*4**2*3**4
lower=Q(1,2**18)*Q(2,3)**12*Q(angle_leading,2**12)*Q(4,9)**4
upper=Q(1,2**18)*Q(2,3)**12*angle_leading*Q(16,9)**4
assert lower==Q(1,688747536)
assert upper==Q(65536,43046721)
assert lower/compound_bound==Q(1,8230118400)
assert upper/Q(5,59049)==Q(65536,3645)
assert Q(1,8230118400)/5==Q(1,41150592000)
assert Q(28)*Q(4)*Q(2)==224  # coefficient of sqrt(theta/2), equivalent to 112 sqrt(2theta)
# Delta >=4g >=8h; throughout h<=1/32, 8h<=1/4<2.
assert 12*(2+Q(1,32))*Q(1,32)==Q(195,256)<1
assert 8*Q(1,32)==Q(1,4)<2
# theta^6 <=theta <= h for theta<=1/8192, hence 8h > theta^6/41150592000.
assert Q(1,8192)<1
# Nominal CP comparison: |H-H0|_max <= (5/3) theta^2.
assert max(Q((n+1)**2,3**n) for n in range(9))==Q(4,3)
# A fixed tilted-axis qubit CP comparison has a rank-four table at O(theta^6).
a=Q(9,34); p=1-a; r2=Q(68,5)
mu={2*k:Q(2**(2*k)+4**(2*k),2) for k in range(5)}
alt={2*k:p*r2**k+(a if k==0 else 0) for k in range(5)}
for k in (0,2,4):assert mu[k]==alt[k]
assert mu[6]-alt[6]==Q(1152,5)
assert (mu[6]-alt[6])/Q(4*factorial(6))==Q(2,25)
cp_remainder=(mu[8]+alt[8])/Q(4*factorial(8))
assert cp_remainder==Q(5669,15750)
cp_bounds=[Q(1,3**n)*(Q(2,25)*(n+1)**6+cp_remainder*(n+1)**8*T) for n in range(9)]
assert max(cp_bounds)<16 and cp_bounds.index(max(cp_bounds))==4

result={
 'arithmetic':'Fraction rationals; SymPy exact polynomial identities; no floating sign samples',
 'requirements':{'python':'>=3.12','sympy':'1.14.0'},
 'source':{'commit':PIN,'path':SOURCE,'sha256':hashlib.sha256(source_bytes).hexdigest()},
 'truncation':'Q[t]/(t^8); t=theta^2',
 's_n_coeff_formula':'[t^k]s_n=1_{k=0}/4+(-1)^k(2^(2k)+4^(2k))(n+1)^(2k)/(8*3^n*(2k)!)',
 'characteristic_coefficients':{str(j):[str(x) for x in cs] for j,cs in char.items()},
 'F_coefficients':[str(x) for x in F],
 'all_4x4_minors':allminors,
 'small_eigenvalue_leading':{'theta^2':str(-A/P),'theta^4':str(-B/A),'theta^6':str(D/B)},
 'schur_identities':schur,
 'determinant_bound':{'lower':str(lower),'upper':str(upper),'power':12},
 'minor_bound':str(minor_bound),'compound_bound':str(compound_bound),
 'sigma_bound':{'lower':'1/8230118400','upper':'65536/3645','power':6},
 'entry_threshold':'theta^6/41150592000',
 'hypothetical_CP_rank4_comparison':{'axis_squared_z':str(a),'rotating_weight':str(p),'angle_squared_factor':str(r2),'sixth_moment_difference':str(mu[6]-alt[6]),'eighth_order_remainder_coefficient':str(cp_remainder),'max_entry_error_bound':'16 theta^6','verified_uniform_factor':str(max(cp_bounds))},
 'newton_trace_coefficients':{str(k):[str(x) for x in cs] for k,cs in tr.items()},
 'CP_uniform_factors':[str(x) for x in cp_bounds],
 'CP_even_moments':{'original':{str(k):str(v) for k,v in mu.items()},'comparison':{str(k):str(v) for k,v in alt.items()}},
 'hankel_error_multiplicities':[1,2,3,4,5,4,3,2,1],
 'remainder_ratios_at_Tmax':{'t4_to_leading':str(r4),'R5_to_leading':str(r5),'both_strictly_below':'1/8'},
 'certified_identities':['moment exponential sum','determinant modulus Laurent identity','all 5 generalized Vandermonde row patterns','characteristic coefficients by principal minors and independent Newton identities','all 25 4x4 minor Taylor polynomials modulo t^8'],
 'assertions':'all passed',
 'proof_boundary':'Analytic inequalities and spectral deductions are ordinary mathematics in C113 sections 37-48; the finite checks are not Lean/kernel verification.'
}
args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n', encoding='utf-8')
print(json.dumps({'status':'PASS','source_sha256':result['source']['sha256'],'leading':result['small_eigenvalue_leading'],'F_t3':str(F[3]),'F_t4':str(F[4]),'remainder_ratios':result['remainder_ratios_at_Tmax'],'certificate':str(args.output.resolve())},indent=2))

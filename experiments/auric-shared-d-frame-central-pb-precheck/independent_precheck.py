"""Independent rational real-coordinate certificate and full-output PB witness checks."""
from itertools import combinations
from hashlib import sha256
import json
import sys
import numpy as np
from flint import fmpq as Q, fmpq_mat as QM, nmod_mat

assert sys.version_info >= (3, 12)
P = 1009
IJ = list(combinations(range(4), 2))
WEDGE = list(combinations(range(15), 2))
ZERO = QM(4, 4)
IDENT = QM([[int(i == j) for j in range(4)] for i in range(4)])

class Gaussian:
    def __init__(self, real, imag=None):
        self.real = real
        self.imag = ZERO if imag is None else imag
    def __matmul__(self, other):
        return Gaussian(self.real * other.real - self.imag * other.imag,
                        self.real * other.imag + self.imag * other.real)
    def star(self):
        return Gaussian(self.real.transpose(), -self.imag.transpose())
    def __eq__(self, other):
        return self.real == other.real and self.imag == other.imag

ONE = Gaussian(IDENT)

def hermitian(diag=(0,0,0,0), real=None, imag=None):
    r = QM([[diag[i] if i == j else 0 for j in range(4)] for i in range(4)])
    s = QM(4,4)
    for (i,j), v in (real or {}).items():
        r[i,j] = r[j,i] = v
    for (i,j), v in (imag or {}).items():
        s[i,j], s[j,i] = v, -v
    return Gaussian(r,s)

BASIS = []
for i in range(3):
    diag = [0]*4
    diag[i], diag[3] = 1,-1
    BASIS.append(hermitian(diag))
for ij in IJ:
    BASIS.extend([hermitian(real={ij:1}), hermitian(imag={ij:1})])

def coordinates(a):
    assert a == a.star()
    tr = sum(a.real[i,i] for i in range(4))
    return [a.real[i,i] - tr/Q(4) for i in range(3)] + [z for i,j in IJ for z in (a.real[i,j], a.imag[i,j])]

def raw_coordinates(a):
    return [a.real[i,i] for i in range(4)] + [z for i,j in IJ for z in (a.real[i,j],a.imag[i,j])]

for k, b in enumerate(BASIS):
    assert coordinates(b) == [Q(int(j == k)) for j in range(15)]
assert coordinates(ONE) == [Q(0)]*15
assert QM([raw_coordinates(b) for b in [ONE]+BASIS]).rank() == 16

def modq(x):
    num, den = int(x.numerator), int(x.denominator)
    assert den % P
    return num * pow(den, -1, P) % P

def modular(a):
    return np.array([[modq(a[i,j]) for j in range(a.ncols())] for i in range(a.nrows())],dtype=np.int64)

seed = 1
denominators = []
def unitary():
    global seed
    u = ONE
    for _ in range(4):
        z = []
        for _ in range(8):
            seed = (1664525*seed + 1013904223) & 0xffffffff
            z.append(seed % 5 - 2)
        r, s = z[:4], z[4:]
        d = sum(x*x for x in z)
        assert 5 <= d <= 26 and d % P
        denominators.append(d)
        hr = QM([[int(i==j) - Q(2*(r[i]*r[j]+s[i]*s[j]),d) for j in range(4)] for i in range(4)])
        hi = QM([[Q(-2*(s[i]*r[j]-r[i]*s[j]),d) for j in range(4)] for i in range(4)])
        h = Gaussian(hr,hi)
        assert h.star() == h and h@h == ONE
        u = u@h
    assert u.star()@u == ONE and u@u.star() == ONE
    # Evaluation and adjoint are distinct maps: compute adjoint before evaluation.
    up = (modular(u.real) + 469*modular(u.imag)) % P
    usp = (modular(u.real.transpose()) - 469*modular(u.imag.transpose())) % P
    assert np.array_equal(up@usp%P,np.eye(4,dtype=np.int64))
    return u

c = hermitian((1,2,3,0))
legs = [hermitian(real={(j,3):1}) if phase == 0 else hermitian(imag={(j,3):1}) for j in (0,1) for phase in (0,1)]
rows, inputs, modular_frames = [], [], []
adjoint_tests = 0
for frame in range(48):
    u = unitary()
    ua = u.star()
    modular_frames.append(((modular(u.real)+469*modular(u.imag))%P, (modular(ua.real)+469*modular(ua.imag))%P))
    forward = QM([coordinates(u@b@ua) for b in BASIS]).transpose()
    backward = QM([coordinates(ua@b@u) for b in BASIS]).transpose()
    assert forward*backward == QM([[int(i==j) for j in range(15)] for i in range(15)])
    assert backward*forward == QM([[int(i==j) for j in range(15)] for i in range(15)])
    assert coordinates(u@c@ua) == list((forward*QM([[a] for a in coordinates(c)])).entries())
    ac = coordinates(u@c@ua)
    output = [[modq(backward[j,k]) for k in range(15)] for j in range(3,15)]
    for t in legs:
        bt = coordinates(u@t@ua)
        # Wedge is formed over Q before reduction, independently of field-i arithmetic.
        w = [modq(ac[i]*bt[j]-ac[j]*bt[i]) for i,j in WEDGE]
        inputs.append(w)
        for functional in output:
            rows.append([wi*dk % P for wi in w for dk in functional])
    # All 105 independent central coefficients survive: output(I)=0 in every frame.
    assert coordinates(ua@ONE@u)[3:] == [Q(0)]*12
    adjoint_tests += 1

m = np.array(rows, dtype=np.int64)
assert m.shape == (2304,1575)
digest = sha256(m.astype('<u2').tobytes()).hexdigest()
assert digest == 'd1b015669a9df51431a8fa8d7e88a977754f0e624760317f28d8b5a900b0afde'
mm = nmod_mat(rows,P)
assert mm.rank() == 1575
rr, rank = mm.transpose().rref()
pivots = [next(j for j in range(rr.ncols()) if rr[k,j]) for k in range(rank)]
determinant = int(nmod_mat([rows[j] for j in pivots],P).det())
assert determinant == 796
assert nmod_mat(inputs,P).rank() == 105
assert all(P % k for k in range(2,32)) and 469**2 % P == P-1
assert len(denominators) == 192
assert sha256(np.array(denominators,dtype='<u2').tobytes()).hexdigest() == '904c8029c79ccf3dfe0b18ee019645e4b02509ba1907c6b6cab3c70ee10c9cc2'
assert sha256(np.array(pivots,dtype='<u2').tobytes()).hexdigest() == '964c5d9faa2479b19edafc7f60b3a9469e5c14da32e7e33b146a6742fbb8f591'

# Test the coefficient tensor convention against a separate direct evaluation.
rng = np.random.default_rng(70623)
coef = rng.integers(0,P,size=(105,15),dtype=np.int64)
flat = coef.ravel()
readings = m@flat % P
basis_mod = [(modular(b.real)+469*modular(b.imag))%P for b in BASIS]
for n,w in enumerate(inputs):
    value = np.array(w,dtype=np.int64)@coef % P
    output_matrix = sum((int(v)*b for v,b in zip(value,basis_mod)),start=np.zeros((4,4),dtype=np.int64))%P
    forward_u,back_u = modular_frames[n//4]
    rotated = back_u@output_matrix%P@forward_u%P
    actual = [z for i,j in IJ for z in (((int(rotated[i,j])+int(rotated[j,i]))*pow(2,-1,P))%P,((int(rotated[i,j])-int(rotated[j,i]))*pow(2*469,-1,P))%P)]
    assert list(readings[n*12:n*12+12]) == actual
print(json.dumps({'independent_rational_reconstruction':'passed', 'rational_unitary_and_adjoint_checks':adjoint_tests,'quotient_basis_rank_with_identity':16,'central_kernel_dimension':105,'matrix_shape':list(m.shape),'matrix_sha256_u16le':digest,'rank_mod_1009':1575,'selected_minor_determinant':determinant,'input_wedge_rank':105,'denominator_count':len(denominators),'denominator_minmax':[min(denominators),max(denominators)],'modular_i_used_for_coordinates':False}),flush=True)

# Derive full-M4 budgets by the real Cauchy--Binet sum of all 2x2 column minors.
import sympy as sp
s = sp.symbols('s',real=True,nonzero=True)

def real_column(v):
    return sp.Matrix([part for z in v for part in (sp.re(sp.expand(z)),sp.im(sp.expand(z)))])

def area2(v,w):
    vr,wr = real_column(v),real_column(w)
    return sp.factor(sum((vr[i]*wr[j]-vr[j]*wr[i])**2 for i,j in combinations(range(8),2)))

J = lambda a,b:(a*b+b*a)/2
identity4=sp.eye(4)
cc=sp.diag(1,2,3,0)
count=0
for a,j in [(1,0),(2,1)]:
    for sigma in (-1,1):
        iota=sigma*sp.I
        for phase in (1,sp.I):
            tt=sp.zeros(4); tt[j,3]=phase; tt[3,j]=sp.conjugate(phase)
            lam= a*(a*a-4*s*s)/(4*s*(a*a+2))
            xx=sp.zeros(4); xx[j,j]=lam
            gg=iota*(cc*tt-tt*cc)/2+s*identity4
            # h=0 follows from this chosen diagonal actor probe.
            assert cc*xx-xx*cc==sp.zeros(4)
            rr=J(cc,xx); kt=J(cc,tt)
            ll=2*J(xx,kt)-2*J(tt,rr+gg)
            v=sp.zeros(4,1); v[j]=phase*s*(1-2*iota/a); v[3]=1
            ell=sp.factor((v.conjugate().T*ll*v)[0])
            readings=[]
            for delta in (-1,1):
                R=xx-delta*kt
                S=tt+delta*(rr+gg)
                readings.append(area2(R*v,S*v))
            if a==1:
                B=(80*s**6+12*s**4+6*s*s+1)/576
                L=-(44*s*s+1)/12
                margin=s*s*(41+962*s*s-40*s**4)/72
            else:
                B=(2*s**6+6*s**4+9*s*s+1)/9
                L=-2*(5*s*s+1)/3
                margin=4*s*s*(1+19*s*s-2*s**4)/9
            assert sp.simplify(ell-L)==0
            assert all(sp.simplify(z-B)==0 for z in readings)
            assert sp.simplify(ell**2-4*B-margin)==0
            assert sp.simplify((v.conjugate().T*v)[0]-(1+s*s*(1+sp.Rational(4,a*a))))==0
            # Large witness is normalized, with x=0. Keep all four output rows.
            z=sp.zeros(4,1); z[j]=phase/sp.sqrt(2); z[3]=1/sp.sqrt(2)
            largeL=-2*J(tt,gg)
            assert sp.simplify((z.conjugate().T*largeL*z)[0]+2*s)==0
            for delta in (-1,1):
                assert sp.simplify(area2(-delta*kt*z,(tt+delta*gg)*z)-sp.Rational(a**4,16))==0
            count+=1
# Algebraic positive bounds for all s in the small intervals, including endpoints.
y=sp.symbols('y',nonnegative=True)
assert sp.expand((41+962*y-40*y*y)-(41+sp.Rational(1919,2)*y)-40*y*(sp.Rational(1,16)-y))==0
assert sp.expand((1+19*y-2*y*y)-(1+17*y)-2*y*(1-y))==0
print(json.dumps({'full_output_real_gram_witness_derivations':'passed','actor_orientation_phase_cases':count,'spectator_output_retained':True,'small_bounds':'a=1: 0<s^2<=1/16; a=2: 0<s^2<=1','large_bound':'abs(s)>a^2/4','strict_positive_budgets':True,'normalization_mass':'any prescribed finite mu>0, scaling by mu/(xi*xi)','zero_mass':'zero positive functional only; not discriminating'}),flush=True)

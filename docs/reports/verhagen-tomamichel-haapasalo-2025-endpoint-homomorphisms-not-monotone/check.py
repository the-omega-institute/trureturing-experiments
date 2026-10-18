#!/usr/bin/env python3
"""Check: the additional homomorphisms of Remark 15 in Verhagen-Tomamichel-Haapasalo (arXiv:2507.07520v3) are not
monotone on S_e.o.: Phi_{alpha,z} (Eq. (84) form) for alpha in {0,1}, z<1, z != 0, in both R_+ and R_+^op, and the
tropical Phi^T_z (Eq. (85) form) for z<0 in TR_+. Part A is exact (sympy): the two Kraus channels, their outputs and
block data. Part B evaluates the homomorphisms on a grid of z. Exit 0 iff all checks pass."""
import sys
import sympy as sp
fails = []
def check(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)
    if not ok: fails.append(name)
s3 = sp.sqrt(3)
A = sp.Matrix([[1, 0], [0, 0]]); B = sp.Matrix([[1, 1], [1, 1]]) / 2
L1 = sp.Matrix([[sp.Rational(1, 2), 0], [0, s3 / 2], [0, 0], [0, 0]])
L2 = sp.Matrix([[0, 0], [0, 0], [s3 / 2, 0], [0, sp.Rational(1, 2)]])
check('flag channel is trace preserving: L1^H L1 + L2^H L2 = I2', sp.simplify(L1.H * L1 + L2.H * L2 - sp.eye(2)) == sp.zeros(2, 2))
Ap = L1 * A * L1.H + L2 * A * L2.H; Bp = L1 * B * L1.H + L2 * B * L2.H
check("A' = diag(1/4, 0, 3/4, 0)", sp.simplify(Ap - sp.diag(sp.Rational(1, 4), 0, sp.Rational(3, 4), 0)) == sp.zeros(4, 4))
Bshown = sp.Matrix([[1, s3, 0, 0], [s3, 3, 0, 0], [0, 0, 3, s3], [0, 0, s3, 1]]) / 8
check("B' = (1/8)[[1,√3,0,0],[√3,3,0,0],[0,0,3,√3],[0,0,√3,1]]", sp.simplify(Bp - Bshown) == sp.zeros(4, 4))
# block data: blocks {0,1} and {2,3}
a = [sp.Matrix([1, 0]), sp.Matrix([1, 0])]; p = [sp.Rational(1, 4), sp.Rational(3, 4)]
b = [sp.Matrix([sp.Rational(1, 2), s3 / 2]), sp.Matrix([s3 / 2, sp.Rational(1, 2)])]; q = [sp.Rational(1, 2), sp.Rational(1, 2)]
blkA = sp.diag(p[0] * a[0] * a[0].H, p[1] * a[1] * a[1].H); blkB = sp.diag(q[0] * b[0] * b[0].H, q[1] * b[1] * b[1].H)
check("A' and B' have the block decomposition (37) with p=(1/4,3/4), q=(1/2,1/2)", sp.simplify(blkA - Ap) == sp.zeros(4, 4) and sp.simplify(blkB - Bp) == sp.zeros(4, 4))
ov = [sp.simplify(abs((a[i].H * b[i])[0]) ** 2) for i in range(2)]
check(f"output squared overlaps are {ov} = (1/4, 3/4), all nonzero (V_e.o.)", ov == [sp.Rational(1, 4), sp.Rational(3, 4)])
check('input pair: unit vectors with squared overlap 1/2 (V_e.o.)', sp.simplify(abs((sp.Matrix([1, 0]).H * sp.Matrix([1, 1]) / sp.sqrt(2))[0]) ** 2) == sp.Rational(1, 2))
E1 = sp.Matrix([[1, 0]]); E2 = sp.Matrix([[0, 1]])
check('erasure channel is trace preserving and sends (A, B) to (1, 1)', sp.simplify(E1.H * E1 + E2.H * E2 - sp.eye(2)) == sp.zeros(2, 2) and (E1 * A * E1.H + E2 * A * E2.H)[0] == 1 and sp.simplify((E1 * B * E1.H + E2 * B * E2.H)[0]) == 1)
def Phi(alpha, z, pp, qq, ovs):
    return sum((pp[i] ** alpha) * (qq[i] ** (1 - alpha)) * ovs[i] ** z for i in range(len(pp)))
def Psi(z, ovs): return max(o ** z for o in ovs)
bad = []
zs = [sp.Rational(k, 20) for k in range(-100, 20) if k != 0]
for alpha in (0, 1):
    for z in zs:
        if alpha == 0:
            fin = Phi(0, z, [1], [1], [sp.Rational(1, 2)]); fout = Phi(0, z, p, q, ov)
        else:
            fin = Phi(1, z, [1], [1], [sp.Rational(1, 2)]); fout = Phi(1, z, q, p, ov)   # swapped pair (B, A) -> (B', A')
        erase = Phi(alpha, z, [1], [1], [1])
        # monotone into R_+ needs Phi(out) <= Phi(in) for every channel; into R_+^op needs Phi(out) >= Phi(in)
        viol_R = (fout > fin) or (erase > fin)
        viol_op = (fout < fin) or (erase < fin)
        if not (viol_R and viol_op): bad.append((alpha, z, float(fin), float(fout), float(erase)))
check(f'Phi_(alpha,z) violates both orders for alpha in (0,1) and all {len(zs)} grid values z in [-5, 1), z != 0', not bad)
if bad: print('INFO first bad', bad[:3])
tbad = [z for z in zs if z < 0 and not (Psi(z, ov) > Psi(z, [sp.Rational(1, 2)]))]
check('tropical Phi^T_z violates TR_+ monotonicity for every grid z < 0', not tbad)
print('INFO z=1/2:', Phi(0, sp.Rational(1, 2), [1], [1], [sp.Rational(1, 2)]), '->', sp.nsimplify(Phi(0, sp.Rational(1, 2), p, q, ov)))
print('INFO z=-1:', Phi(0, -1, [1], [1], [sp.Rational(1, 2)]), '->', Phi(0, -1, p, q, ov), '; tropical', Psi(-1, [sp.Rational(1, 2)]), '->', Psi(-1, ov))
print('RESULT', 'FAIL' if fails else 'PASS', len(fails))
sys.exit(1 if fails else 0)

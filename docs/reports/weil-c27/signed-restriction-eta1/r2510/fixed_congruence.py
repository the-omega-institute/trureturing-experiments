"""Finite midpoint-selected DYADIC BASIS; signs require outward full112 replay.

The basis is exactly invertible by triangularity, not by approximate accuracy.
No midpoint or approximate factor sign is used as mathematical evidence.
"""
from decimal import Decimal, localcontext, ROUND_FLOOR
import hashlib
import json

REDUCTION_SHA256='b46df8ada7963dc098e1cf92809c3081ed31017b261f3b189009f156fc6ee4e4'
BITS=4096


def _point(value, exponent, bits):
    numerator,denominator=value.as_integer_ratio()
    shift=bits+exponent
    if shift>=0:numerator<<=shift
    else:denominator<<=-shift
    return numerator//denominator


def propose_basis(raw):
    assert hashlib.sha256(raw).hexdigest()==REDUCTION_SHA256
    packet=json.loads(raw)
    assert packet['arithmetic_bits']==1024 and packet['reduced_dimension']==112
    with localcontext() as ctx:
        ctx.prec=300
        n=112
        a=[[None]*n for _ in range(n)]
        divisor=Decimal(2)**1025
        for i,row in enumerate(packet['reduced112_upper']):
            for off,(lo,hi) in enumerate(row):
                j=i+off
                a[i][j]=a[j][i]=(Decimal(int(lo))+Decimal(int(hi)))/divisor
        row_exp=[0]*n
        ln2=Decimal(2).ln()
        for _ in range(8):
            steps=[]
            for row in a:
                largest=max(abs(x) for x in row)
                log2=int((largest.ln()/ln2).to_integral_value(rounding=ROUND_FLOOR)) if largest else 0
                steps.append(-(log2//2))
            for i in range(n):
                row_exp[i]+=steps[i]
                for j in range(i,n):
                    a[i][j]=a[j][i]=a[i][j]*(Decimal(2)**(steps[i]+steps[j]))
        u=[[Decimal(int(i==j)) for j in range(n)] for i in range(n)]
        col_exp=[None]*n
        active=list(range(n))
        blocks=[]
        while active:
            m=len(active)
            di=max(range(m),key=lambda i:abs(a[i][i]))
            diag=abs(a[di][di]);off=Decimal(0);pair=None
            for i in range(m):
                for j in range(i+1,m):
                    v=abs(a[i][j])
                    if v>off:off=v;pair=[i,j]
            assert max(diag,off)>0, 'Approximate zero suffix is not a singularity certificate.'
            pivot=[di] if pair is None or diag>=off/2 else pair
            original=[active[i] for i in pivot]
            blocks.append(original)
            largest=max(abs(a[i][j]) for i in pivot for j in pivot)
            exponent=-(int((largest.ln()/ln2).to_integral_value(rounding=ROUND_FLOOR))//2)
            for label in original:col_exp[label]=exponent
            if len(pivot)==1:inverse=[[1/a[di][di]]]
            else:
                i,j=pivot;aa,bb,dd=a[i][i],a[i][j],a[j][j]
                det=aa*dd-bb*bb;assert det!=0
                inverse=[[dd/det,-bb/det],[-bb/det,aa/det]]
            rest=[i for i in range(m) if i not in pivot]
            cross=[[a[i][k] for k in pivot] for i in rest]
            mult=[[sum((row[k]*inverse[k][j] for k in range(len(pivot))),Decimal(0)) for j in range(len(pivot))] for row in cross]
            # Fixed approximate shears, cached pivot columns; no sign is certified.
            pivot_columns=[[u[i][label] for label in original] for i in range(n)]
            for j,old_j in enumerate(rest):
                label=active[old_j]
                for i in range(n):
                    u[i][label]-=sum((pivot_columns[i][k]*mult[j][k] for k in range(len(pivot))),Decimal(0))
            updated=[[None]*len(rest) for _ in rest]
            for i,old_i in enumerate(rest):
                for j in range(i,len(rest)):
                    value=a[old_i][rest[j]]-sum((mult[i][k]*cross[j][k] for k in range(len(pivot))),Decimal(0))
                    updated[i][j]=updated[j][i]=value
            active=[active[i] for i in rest];a=updated
        order=[i for block in blocks for i in block]
        assert sorted(order)==list(range(n))
        basis=[]
        for i,ri in enumerate(order):
            assert u[ri][ri]==1 and all(u[ri][order[j]]==0 for j in range(i))
            row=[]
            for j in range(i,n):
                cj=order[j]
                row.append(str(_point(u[ri][cj],row_exp[ri]+col_exp[cj],BITS)))
            assert int(row[0])>0
            row_exponent=row_exp[ri]+col_exp[ri]
            assert BITS+row_exponent>=0 and int(row[0])==1<<(BITS+row_exponent)
            basis.append(row)
    return dict(schema='R2510_FIXED_INVERTIBLE_DYADIC_BASIS_V1',
        purpose='inertia_research_only',reduction_result_sha256=REDUCTION_SHA256,
        original_arithmetic_bits=1024,basis_grid_bits=BITS,dimension=112,
        permutation=order,original_label_blocks=blocks,
        basis_upper_point_integers=basis,
        basis_encoding='upper triangular in permuted original coordinates; exact integer/2^basis_grid_bits',
        invertibility='strictly positive exact power-of-two diagonal; determinant is their product',
        midpoint_precision_decimal_digits=300,midpoints_select_basis_only=True,
        expected_inertia=None)


def enclose_congruence(core, base, matrix, upper_point_integers, bits):
    """Parameter dimension permits independent small rational product checks."""
    n=len(matrix)
    core.require(n>0 and len(upper_point_integers)==n and all(len(row)==n for row in matrix),
                 'shape_mismatch','congruence dimension')
    zero=core.Interval(0,0,bits)
    basis=[[zero for _ in range(n)] for _ in range(n)]
    for i,row in enumerate(upper_point_integers):
        core.require(len(row)==n-i,'shape_mismatch','upper dyadic basis')
        for off,raw in enumerate(row):
            core.require(type(raw) is str and str(int(raw))==raw,'invalid_basis','canonical exact integer')
            point=int(raw);basis[i][i+off]=core.Interval(point,point,bits)
        point=int(row[0])
        core.require(point>0 and point&(point-1)==0,'singular_basis','positive power-of-two diagonal required')
    lifted=[[base._lift(core,x,bits) for x in row] for row in matrix]
    core.require(all(lifted[i][j]==lifted[j][i] for i in range(n) for j in range(n)),
                 'invalid_symmetric_state','same symmetric realization')
    right=[[zero for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            value=zero
            for k in range(j+1):
                if not basis[k][j].is_zero():value=value+lifted[i][k]*basis[k][j]
            right[i][j]=value
    out=[[zero for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(i,n):
            value=zero
            for k in range(i+1):
                if not basis[k][i].is_zero():value=value+basis[k][i]*right[k][j]
            out[i][j]=out[j][i]=value
    return out

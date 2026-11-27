"""在五种荷的双钟环面上分别认证三条边界射线首站的精确盲点。"""
from pathlib import Path
from runpy import run_path
from flint import acb, acb_mat, arb, arb_mat, ctx


CASES = ((1, 2, ".43052262", "5.82209977"), (-4, 1, "5.60000358", ".87833265"),
         (2, 2, ".55974120", "5.64213905"), (-3, 1, ".63870109", "5.50779883"),
         (6, 2, "5.54725202", "1.03239116"))


def verify(root):
    graph = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    pump = run_path(str(root/"ch39-native-light-projector-pump.py"))
    ctx.prec=160
    frames=list(pump["rotation_bases"](graph))
    sites=graph["geo"]["sites"]
    boundary=sorted(graph["boundary_sites"]())
    assert boundary==[6,7,8]
    assert {sites.index((y,z,x)) for x,y,z in (sites[i] for i in boundary)}==set(boundary)
    for q,character,first,second in CASES:
        base,links=pump["clock_matrices"](graph,q)
        adj=pump["adjoint"]
        frame=frames[character]
        fixed=adj(frame)*base*frame
        blocks=[adj(frame)*m*frame for m in links]
        n=fixed.nrows()
        eye=acb_mat([[int(i==j) for j in range(n)] for i in range(n)])
        target=[j for j in range(n) if abs(frame[boundary[0],j])>arb(".5")]
        assert target==[2]
        inactive=target[0]

        def matrix_angles(a,b):
            terms=[block*acb(0,angle).exp() for block,angle in zip(blocks,(a,b),strict=True)]
            return fixed+sum((t+adj(t) for t in terms),acb_mat(n,n)), [acb(0,1)*(t-adj(t)) for t in terms]

        origin=[arb(first),arb(second)]
        initial,_=matrix_angles(*origin)
        v,gap=pump["first_line"](initial)
        anchor=max(range(n),key=lambda k:abs(v[k]).lower())
        assert anchor!=inactive and gap>arb(1)
        rest=[k for k in range(n) if k not in (inactive,anchor)]
        normalized=[z/v[anchor] for z in v]
        ev=(adj(acb_mat([[z] for z in v]))*initial*acb_mat([[z] for z in v]))[0,0].real
        center=origin+[ev.mid()]+[a.mid() for k in rest for a in (normalized[k].real,normalized[k].imag)]
        dim=len(center)
        assert dim==2*n-1

        def equations(values, include_soft_error):
            mat, deriv=matrix_angles(values[0],values[1])
            # 镜面对称的三条辐边同步改权，19站扰动范数 <= 6e-8。
            # 两个钟角独立变化，软辐边参数固定；区间覆盖每个允许取值。
            if include_soft_error:
                entry=acb(arb(0,arb("6e-8")),arb(0,arb("6e-8")))
                mat+=acb_mat([[entry for _ in range(n)] for _ in range(n)])
            shifted=mat-eye*values[2]
            vector=acb_mat(n,1)
            vector[anchor,0]=1
            for idx,k in enumerate(rest):
                vector[k,0]=acb(values[3+2*idx],values[4+2*idx])
            residual=shifted*vector
            derivatives=[op*vector for op in deriv]+[-vector]
            for k in rest:
                line=acb_mat([[shifted[i,k]] for i in range(n)])
                derivatives.extend((line,acb(0,1)*line))
            def unpack(column):
                return [component for i in range(n) for component in
                        ((column[i,0].real,) if i==anchor else (column[i,0].real,column[i,0].imag))]
            return arb_mat([[z] for z in unpack(residual)]),arb_mat([unpack(d) for d in derivatives]).transpose()

        radius=arb(".0001")
        box=[x+arb(0,radius.upper()) for x in center]
        _,jac=equations(center,False)
        pre=jac.mid().inv()
        assert not pre.det().contains(0)
        offset,_=equations(center,True)
        _,jacbox=equations(box,True)
        I=arb_mat([[int(i==j) for j in range(dim)] for i in range(dim)])
        defect=I-pre*jacbox
        displacement=pre*offset
        ratio=max(abs(displacement[i,0]).upper()/radius+
                  sum(abs(x).upper() for x in defect.tolist()[i]) for i in range(dim))
        assert ratio<1,(q,ratio)
        # 其余虚分量由 Hermitian 二次型恒等式恢复。全盒根仍在字符最轻本征值的隔离隙中。
        shift=arb(2).sqrt()*2*radius+arb("6e-8")
        assert box[2].upper()<ev.lower()+gap-shift
        print("荷、字符、固定坐标、钟角盒半宽、全盒自映比",q,character,anchor,
              radius.str(4),ratio.str(11),flush=True)


if __name__=="__main__":
    import sys
    root=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
    verify(root)

"""两个三项实铅笔的正谱函数同时满足六质量窗和最小 CKM 混合角。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb


def verify():
    replay = Path(__file__).resolve().parent
    model,up,down = run_path(str(replay/"ch39-two-zero-orbit-condition.py"))["calibration"]()
    frame,h = run_path(str(replay/"ch39-shell-pencil-six-mass-window.py"))["shell_frame"](model)
    ratio = h[0,1].real/h[1,2].real
    adj = model["adj"]
    cases = (
        (("1.029923870182003 .0006721828059758628 .00003756287337657452",
          ".3257568433117072 -.004388276603562354 -.000005656289276394763",
          ".7365313191751384 .0010009269093679 -.00005002411170360671"),
         "19.596131826395034","2.7723157267597989"),
        ((".9788401815713665 .001182835895828091 .00009695629380393149",
          ".5161383481948316 -.002681451490776185 .00002409262405821096",
          ".696284559591118 .0003248566079149002 -.0001541608842793307"),
         "15.835644452588557","2.2403097993123082"),
        (("-1.0299955061154296 .0006739588194340176 -.00003734351146289303",
          "-.325480640747518 -.004401646063320289 .000005612213663910973",
          "-.7365532522102924 .0010026105016756317 .000049741117807932745"),
         "19.604138072394083","2.774685099285122"),
    )
    for index,(rows,su,sd) in enumerate(cases):
        factor = acb_mat([[arb(v) for v in row.split()] for row in rows])
        assert abs(factor.det())>0
        levels,basis = model["eig"](factor*adj(factor))
        assert 0 < levels[0] < levels[1] < levels[2]
        for i in (1,2):
            y,z = [(basis[j,i]/basis[0,i]).real for j in (1,2)]
            first = y*(1-ratio*z)
            second = z*(ratio+z)-y*y
            eigenvalue = ratio+z-ratio*y*y
            assert abs((ratio*second-eigenvalue)*(-eigenvalue)-first**2)>0
        matrices = [arb(scale)*m*frame*factor for scale,m in ((su,up),(sd,down))]
        (u,U),(d,D) = [model["eig"](adj(m)*m) for m in matrices]
        values = [v.sqrt() for v in u+d]+model["mixing_angles"](U,D)
        errors = [abs(v-m)/t for v,m,t in zip(values,model["ref"],model["tol"],strict=True)]
        matched = (0,1,2,3,4,5,8) if index<2 else (0,1,2,3,4,5,6,8)
        held_out = [i for i in range(10) if i not in matched]
        bound = arb(".0014") if index<2 else arb(".096")
        assert all(errors[i] < bound for i in matched)
        assert all(errors[i] > 1 for i in held_out)
        print("入窗项数、最大误差占半宽",len(matched),max(errors[i].upper() for i in matched).str(10))
        print("四项混合读数",[v.str(11) for v in values[6:]])
        print("留出项误差占半宽",[errors[i].str(8) for i in held_out])


if __name__=="__main__":
    verify()

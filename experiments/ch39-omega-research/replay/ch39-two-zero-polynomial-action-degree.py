"""用真实零点与原相位的单位辐边极限检验共用的低阶多项式谱作用。"""
from pathlib import Path
from runpy import run_path
from flint import acb, arb, acb_mat


def verify():
    replay = Path(__file__).resolve().parent
    selector = run_path(str(replay/"ch39-light-spectrum-zero-selector.py"))
    light = selector["light_spectra"](replay)
    moments = selector["chebyshev_moments"]
    radius = selector["interior_radius"]
    centers = {-4:2345666, 2:964890}
    shell = run_path(str(replay/"ch39-three-ray-chiral-boundary.py"))
    pump = run_path(str(replay/"ch39-native-light-projector-pump.py"))
    radial = acb_mat(19,19)
    pole = shell["soft"]["pole"]
    for i in shell["soft"]["star"]["orbit"]:
        radial[pole,pole] += 1
        radial[i,i] += 1
        radial[pole,i] -= 1
        radial[i,pole] -= 1

    def flat_energies(q):
        h = acb.zeta_zero(centers[q]).imag
        delta = arb(23)/25*arb(2)**-20/(h*h+shell["soft"]["point"]["study"]["offset"])
        flat = shell["laplacian"](q,h)+delta*radial
        result = []
        for block in pump["rotation_bases"](shell):
            spec = (pump["adjoint"](block)*flat*block).eig(algorithm="rump")
            assert all(z.is_finite() and z.imag.contains(0) for z in spec)
            levels = sorted((z.real for z in spec),key=lambda a:float(a.mid()))
            assert levels[0]>0 and all(levels[0]<v for v in levels[1:])
            result.append(levels[0])
        return result,delta

    finite = ((-4,31681),(-4,1162532),(2,239010),(2,500730),(2,2634783),(2,23551))
    closure = ((-4,None),(-4,783632),(-4,1162532),(-4,8551),(-4,14801),(-4,24501),
               (2,None),(2,2028056),(2,8751),(2,11751))
    for points in (finite,closure):
        degree = len(points)-1
        vectors = []
        for q,n in points:
            if n is None:
                candidate,scale = flat_energies(q)
                assert scale>0
            else:
                candidate,scale = light(q,n),arb(1)
            base = moments(light(q,centers[q]),degree)
            vectors.append([(b-a)/scale for a,b in zip(base,moments(candidate,degree),strict=True)])
        r = radius(vectors)
        assert r > arb("5.37486226e-7" if degree == 5 else "1.2650e-11")
        print("共同谱多项式次数、按所列正比例缩放后的矩位移内半径",degree,r.str(16))


if __name__=="__main__":
    verify()

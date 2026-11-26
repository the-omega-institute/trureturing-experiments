"""在原十九站壳的实际零点间构造长面主值整数及其三角细分。"""
from functools import cache
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx


TRIPLE = (22, 114, 70)
QUAD = (2511, 4648, 5413, 3385)
DISTINCT = (
    2511, 4648, 5413, 3385, 3730, 789, 4886, 817, 2885, 4529, 5209, 5252, 3536, 726, 4088, 289,
    1439, 3934, 1933, 163, 3651, 635, 3928, 3003, 5945, 346, 3501, 19, 2401, 1549, 1511, 4646,
    4520, 2261, 4125, 601, 846, 4931, 2553, 2891, 4802, 948, 1272, 2333, 5042, 2817, 5331, 2927,
    5328, 3388, 5992, 976, 2072, 2080, 1340, 1970, 1612, 165, 818, 4446, 2961, 238, 2739, 5745,
    1509, 1723, 4087, 2111, 2476, 4811, 2369, 881, 1754, 5013, 2404, 753, 876, 214, 1441, 116,
    396, 2519, 3774, 2151, 4359, 487, 4045, 3271, 3457, 5954, 5127, 3615, 4840, 4014, 4924, 3932,
    4482, 5216, 1240, 4929, 4884, 2781, 3460, 3463, 4558, 1479, 2588, 4564, 5367, 5172, 2329, 5951,
    1037, 1142, 262, 2185, 2586, 4849, 3423, 3422, 1965, 1795, 3848, 691, 5289, 2370, 1008, 4327,
    1999, 1309, 848, 4602, 3115, 5664, 4843, 2442, 5699, 5459, 4805, 787, 3306, 1013, 4366, 3810,
    657, 402, 2851, 1275, 51, 1378, 5045, 1339, 4283, 1346, 3004, 2220, 1370, 2118, 1544, 2774,
    2184, 3351, 4285, 5412, 3499, 1620, 5452, 2039, 2546, 372, 1899, 4684, 5004, 3622, 5908, 1860,
    975, 4610, 5254, 5617, 1825, 5748, 2255, 847, 5577, 1587, 2811, 2587, 17, 5870, 3891, 1108,
    5989, 2705, 5374, 4408, 4203, 5420, 455, 1171, 2253, 4292, 316, 427, 3767, 606, 3575, 4129,
    661, 4491, 2078, 1755, 5906, 3159, 2624, 4484, 4401, 1075, 371, 5291, 508, 603, 2964, 4967,
    4006, 5703, 1208, 1547, 4758, 5498, 4761, 4727, 5081, 4969, 18, 5495, 3187, 3892, 2892, 1408,
    311, 3702, 5248, 3927, 3265, 1939, 2629, 5700, 4563, 5051, 5658, 662, 5532, 3581, 5166, 5946,
    5900, 1414, 4447, 2005, 5822, 5133, 164, 1826, 5003, 4131, 2147, 2962, 2179, 5915, 5987, 482,
    3846, 3508, 2001, 2328, 3693, 1935, 4565, 5330, 5778, 913, 1377, 2855, 2626, 2665, 4168, 1865,
    314, 2479, 5336, 5784, 1006, 1762, 5292, 37, 137, 759, 5824, 3542, 2844, 2187, 5746, 1243,
    5450, 2192, 4205, 5661, 1306, 3855, 4963, 3970, 2327, 2746, 114, 2482, 4165, 1146, 912, 3193,
    2813, 139, 3579, 4890, 1068, 5909, 55, 2407, 2661, 4887, 4328, 5990, 4278, 1241, 2589, 3657,
    3033, 399, 784, 5213, 5572, 4567, 1043, 3652, 4080, 1410, 3811, 5912, 3609, 3971, 1204, 4765,
    2657, 3544, 3076, 2186, 1338, 668, 4527, 3773, 5533, 3196, 634, 4440, 32, 263, 3965, 568,
    3416, 5376, 5253, 3267, 1202, 5381, 4083, 3039, 1170, 4174, 4050, 3426, 2997, 5177, 3654, 5126,
    1931, 5008, 3543, 423, 1753, 5211, 1137, 3152, 1135, 5298, 1794, 3505, 4845, 4486, 1827, 5010,
    2441, 4572, 2968, 3044, 3772, 3201, 541, 4923, 4999, 5337, 5089, 4209, 723, 2929, 7, 1759,
    1270, 4686, 5050, 3194, 1787, 3006, 5083, 1510, 3614, 1862, 2780, 2224, 1374, 3232, 1823, 3537,
    3494, 3618, 3304, 5456, 453, 2857, 2403, 4808, 1301, 2405, 2364, 1038, 3229, 2741, 2699, 5417,
    4085, 982, 5619, 3121, 782, 1690, 1577, 1932, 159, 2821, 3849, 1307, 3224, 3429, 2664, 186,
    2219, 1653, 4881, 941,
)


def verify(root):
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    light = run_path(str(root / "ch39-native-light-projector-pump.py"))
    ctx.prec = 192
    basis = list(light["rotation_bases"](shell))[1]
    adj = light["adjoint"](basis)

    @cache
    def line(n):
        full = shell["laplacian"](-3, acb.zeta_zero(n).imag)
        v, _ = light["first_line"](adj * full * basis)
        return v

    @cache
    def link(i,j):
        a, b = line(i), line(j)
        z = sum(v*u.conjugate() for v,u in zip(a,b))
        assert abs(z)**2 > arb(49)/81
        return z/abs(z)

    @cache
    def triangle(i,j,k):
        phase = (link(i,j)*link(j,k)*link(k,i)).arg()
        assert abs(phase) < arb(".38")
        return phase

    def sphere(boundary, target, front, back, name):
        edges = list(zip(boundary, boundary[1:]+boundary[:1]))
        fine = sum(triangle(front,i,j) for i,j in edges)
        opposite = sum(triangle(back,j,i) for i,j in edges)
        factors = [link(j,i) for i,j in edges]
        while len(factors) > 1:
            if len(factors) % 2:
                factors.append(acb(1))
            factors = [a*b for a,b in zip(factors[::2], factors[1::2])]
        coarse = factors[0]
        rounded = (fine+coarse.arg())/(2*arb.pi())
        refined = (fine+opposite)/(2*arb.pi())
        assert rounded.unique_fmpz() == target and refined.unique_fmpz() == 0
        assert (2*target-1)*arb.pi() < fine < (2*target+1)*arb.pi()
        print(name, "边数、正向三角和、反向长面主辐角、粗细两个整数",
              len(boundary), fine.str(16), coarse.arg().str(16), target, 0)

    # 在三标签例中两极允许使用同一个背景，两个网格顶点仍不同。
    print("三枚实际零点及三角相位", TRIPLE, triangle(*TRIPLE).str(17))
    for repetitions,target in ((34,1),(168,3)):
        sphere(TRIPLE*repetitions, target, TRIPLE[0], TRIPLE[1], "三标签")
    sphere(QUAD*22+QUAD[:3], 1, QUAD[0], QUAD[1], "四标签")
    sphere(QUAD*111, 3, QUAD[0], QUAD[1], "四标签")
    # 同一四标签中任意顺序及回返均计入。把每条重叠边的实相位向上取整。
    precision = 10**12
    upper = [[int((precision*triangle(QUAD[0],i,j)).upper().ceil().unique_fmpz())
              for j in QUAD] for i in QUAD]
    distances = upper
    for length in range(1,444):
        cap = arb(max(distances[i][i] for i in range(4)))/precision
        if length < 91:
            assert cap < arb.pi()
        assert cap < 5*arb.pi()
        distances = [[max(distances[i][j]+upper[j][k] for j in range(4))
                      for k in range(4)] for i in range(4)]
    print("指定四枚的两个最短闭边数",91,444)
    assert len(DISTINCT)==len(set(DISTINCT))==468 and {1,2}.isdisjoint(DISTINCT)
    independent = sorted((abs(line(n)[0])**2 for n in (*DISTINCT,1,2)),key=float)
    assert all(a < b for a,b in zip(independent,independent[1:]))
    sphere(DISTINCT[:91],1,1,2,"各顶点零点和投影全异")
    sphere(DISTINCT,3,1,2,"各顶点零点和投影全异")


if __name__ == "__main__":
    root = Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
    del sys.argv[1:]
    verify(root)

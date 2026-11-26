"""十二个站点可同时给出共同标尺的顶、底、粲三项汤川值。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


SITES = [0,1,4,5,7,10,11,12,14,15,17,18]
COEFFICIENTS = [(15689142819,0),(-16744486632,-48489807545),(-9337983568,-23161682277),
                (12014764300,19743072697),(-100868897418,168495516948),(-170116290351,-267503732301),
                (-58704523025,220739875977),(314116292134,411328167501),(38363993806,-417998499912),
                (-432102372162,75751909592),(147745820897,-329895207400),(115006830452,156004850322)]


def verify():
    root = Path(__file__).parent
    flavor = run_path(str(root/"ch39-minimum-flavor-support.py"))
    singulars = run_path(str(root/"ch39-ten-site-charm-exclusion.py"))["singulars"]
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 512
    maps = flavor["flavor_maps"](shell,11431)
    up,down = [acb_mat([[m[i,j] for j in SITES] for i in range(9)]) for m in maps]
    adjoint = flavor["adjoint"]
    identity = acb_mat([[int(i==j) for j in range(12)] for i in range(12)])
    shape = acb_mat([[acb(a,b)/10**12] for a,b in COEFFICIENTS])
    unit = acb_mat([[int(i==0)] for i in range(12)])
    # 两个核中的方向分别改变上型与下型。
    def kernel_direction(matrix):
        right = adjoint(matrix)
        projection = identity-right*(matrix*right).inv()*matrix
        vector = projection*unit
        return vector/sum(abs(row[0])**2 for row in vector.tolist()).sqrt()
    step = kernel_direction(down)*arb(".023")
    bottom_direction = kernel_direction(up)
    assert singulars(down*bottom_direction)[2] > arb(".00001")
    # 两个极小奇异值无需彼此分离。二阶余子式给 s1*s2 的上界。
    initial = up*shape
    matrix = acb_mat([[initial[3*i+j,0] for j in range(3)] for i in range(3)])
    total = sum(abs(z)**2 for row in matrix.tolist() for z in row)
    pairs = list(combinations(range(3),2))
    exterior = sum(abs(matrix[i,k]*matrix[j,l]-matrix[i,l]*matrix[j,k])**2
                   for i,j in pairs for k,l in pairs).sqrt()
    second_upper = exterior/(total/3).sqrt()
    first_lower = (total-2*second_upper**2).sqrt()
    after = singulars(up*(shape+step))
    bottom = singulars(down*shape)[2]
    change = singulars(up*step)[2]
    assert second_upper/first_lower < arb("1e-8")
    assert after[1]/after[2] > arb(".04")
    ratio_floor = (first_lower-change)/bottom
    assert ratio_floor > 70
    print("插值两端的粲顶比、全段顶底比下界",
          (second_upper/first_lower).str(8),(after[1]/after[2]).str(12),ratio_floor.str(12))


if __name__ == "__main__":
    verify()

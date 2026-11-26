"""平行双基本中介在酉背景上的弱矢量质量谱。"""
from sympy import Matrix, cancel, diag, eye, expand, factor, symbols


def verify():
    left, right, x, y, t = symbols("g_L g_R x y t", real=True)
    z = symbols("z")
    cosine = (1 - t**2) / (1 + t**2)
    sine = 2*t / (1 + t**2)
    rotation = Matrix([[cosine, -sine, 0], [sine, cosine, 0], [0, 0, 1]])
    assert (rotation.T * rotation - eye(3)).applyfunc(cancel).is_zero_matrix

    def twice_mass(orientation):
        blocks = [[left**2 * eye(3), -left*right * orientation],
                  [-left*right * orientation.T, right**2 * eye(3)]]
        return Matrix.vstack(*(Matrix.hstack(*row) for row in blocks))

    for orientation in (rotation, diag(-1, -1, 1)):
        c = orientation[0, 0]
        matrix = x * twice_mass(eye(3)) + y * twice_mass(orientation)
        # 先按 x,y 展开，避免把有理圆参数的分母混入系数提取。
        mixed = expand((matrix * matrix).trace() / 4).coeff(x).coeff(y)
        expected = 3*(left**4 + right**4)/2 + left**2*right**2*(1 + 2*c)
        assert cancel(mixed - expected) == 0
        # 块行列式的对角块为标量。六阶多项式只需一个三阶行列式。
        shift = (z - (x+y)*left**2) * (z - (x+y)*right**2)
        coupling = matrix[:3, 3:]
        reduced = shift*eye(3) - coupling*coupling.T
        axial = z * (z - (x+y)*(left**2 + right**2))
        splitting = 2*left**2*right**2*x*y*(1-c)
        target = diag(axial + splitting, axial + splitting, axial)
        assert (reduced-target).applyfunc(cancel).is_zero_matrix
        numerical = {left: 2, right: 3, x: 5, y: 7, t: 2}
        direct = matrix.subs(numerical).charpoly(z).as_expr()
        assert factor(direct - (axial*(axial+splitting)**2).subs(numerical)) == 0

    print("平行中介的混合迹、六阶特征多项式、有理圆及半转端点均通过。")


if __name__ == "__main__":
    verify()

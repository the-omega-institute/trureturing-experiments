# 实际正零点上的十五站满秩与近盲频率

在十九站壳上，固定任意十五站，读取五荷 $q=1,-4,2,-3,6$ 各自独立三维轻态的总复振幅。使用共同站点与相位参考，输入各循环特征的最低线均取单位向量，输出取站点坐标的欧氏范数。以 $s_S(x,\delta)$ 表示布局 $S$ 的最小奇异值。两只物理钟在正黎曼零点纵坐标 $\gamma$ 的角度为 $x_\gamma=(a\gamma,b\gamma)\bmod 2\pi$，辐边软化按模型取 $\delta_\gamma=O(\gamma^{-2})$。以下密度以 $N(T)=\#\{0<\gamma\le T\}$ 为分母，按重数计入所有非平凡零点。

全部 $3876$ 种布局在同一组自然密度为一的实际零点上满秩。另一方面，每个固定布局都会在正比例的零点上接近失秩。若 $\varepsilon(T)\downarrow0$ 且 $\varepsilon(T)\sqrt{\log T}\to\infty$，甚至要求其余布局在同一零点也满秩，仍有
$$
\liminf_{T\to\infty}
\frac{\#\{0<\gamma\le T : s_S(x_\gamma,\delta_\gamma)<\varepsilon(T),\quad
                    \min_U s_U(x_\gamma,\delta_\gamma)>0\}}
     {N(T)\,\varepsilon(T)^2}
\ \ge\ \frac{1}{264\pi^2}.
$$
各布局的近盲零点可不同。是否有实际零点落在某个精确盲点，及各布局在整个环面上的近盲比例首项，仍待判定。[连续钟角上的十五站盲点](ch39-five-charge-fifteen-masks.md)

## 一个能同时读出全部布局的钟角

取编号 $2\,306\,888\,100\,085$ 的正零点。固定同一对钟角，分别用零软化钟边矩阵和实际软化的物理拉普拉斯算子构造十五列。两种情况下，逐一检查 $\binom{19}{15}=3876$ 个方阵 $M$，均得到 $s_S>3.11736\times10^{-6}$。具体地，取中点逆的中点为 $L$，以复球算术给出 $(1-\|I-LM\|_F)/\|L\|_F$ 的下界。正的分子及 $\|LMv\|\le\|L\|_F\|Mv\|$ 保证该量不超过 $s_S$。行枚举直接覆盖每一种布局。[全部布局的复算程序](ch39-five-charge-zero-density.py) · [钉版物理拉普拉斯算子](https://github.com/the-omega-institute/trureturing-experiments/blob/cb58c8e64ecea67184e0c1d68a690be098a7a384/experiments/ch39-omega-research/replay/ch39-three-ray-chiral-boundary.py)

记 $R_S$ 为取行映射，$P_{qk}$ 为嵌入十九站的荷 $q$、循环特征 $k$ 的最低线正交投影。读出矩阵的行 Gram 矩阵为 $G_S=R_S(\sum_{q,k}P_{qk})R_S^*$。每个特征内的谱隙全域大于一，局部围住最低本征值的预解式积分表明 $P_{qk}$ 实解析。这给出与本征向量相位选择无关的实解析函数 $F(x)=\prod_S\det G_S(x,0)$。上面的零软化共同锚点使 $F$ 非恒零，故其零集在连通双钟环面上的 Haar 面积为零。[各特征的钉版谱隙](https://github.com/the-omega-institute/trureturing-experiments/blob/cb58c8e64ecea67184e0c1d68a690be098a7a384/experiments/ch39-omega-research/replay/ch39-native-light-projector-pump.py) · [实解析函数零集的零测性](https://doi.org/10.1134/S0001434620030189)

## 从钟角转到零点纵坐标

令 $\theta_d=\arccos(1-d/1734)$，$a=\theta_{724}-\theta_{244}$，$b=\theta_{804}-\theta_{244}$。三个虚二次根的平方类各有独有的素因子，故可独立共轭；$2\cos\theta_{724}$ 与 $2\cos\theta_{804}$ 都是非整数有理数。于是对 $m\in\mathbb Z^2\setminus\{0\}$，$e^{iu}\ne1$，其中 $u=m_1a+m_2b$。否则独立共轭就会让一个非零幂的相位成为单位根。对每个固定 $u$，Ford–Meng–Zaharescu 引理 1、2 与 $N(T)\asymp T\log T$ 给出 $e^{iu\gamma}$ 的归一化平均趋零，Weyl 判据从而给出两钟的无条件等分布。[零点和的原公式](https://arxiv.org/abs/1511.06814v2)

因 $\delta_\gamma\to0$，在紧环面上，实际读出的奇异值与零软化值之差一致趋零。对任意 $r(T)\downarrow0$，有某个实际布局满足 $s_S\le r(T)$ 的上密度，受每个固定 $\eta>0$ 对应的闭集 $\{\min_S s_S(x,0)\le\eta\}$ 的 Haar 测度控制。令 $\eta\downarrow0$ 就得到零，也包含了精确失秩的情形。对固定布局，记 $M_S=\max_x s_S(x,0)$。若 $0<t<M_S$，在一个 $s_S>t$ 的点上 $\det(G_S-t^2I)\ne0$，故水平集 $s_S=t$ 也是零测集。实际零点中 $s_S<t$ 的极限比例等于平坦环面的相应面积。连续性、已认证盲点和最大点保证此比例严格介于零和一。

## 近盲区域有多大

每类钟边从六个不同的源点连向另三个点，每个目标两边。单个角导数的非对角块具有正交的目标行，各行范数为 $\sqrt2|q|$；源和目标不交，整个自伴导数的算子范数也是 $\sqrt2|q|$。软化导数是三边星形图拉普拉斯的负数，谱在 $[-4,0]$ 中。加上 $2I$ 后范数至多为二，此标量平移不改变本征线。沿任一路径给单位最低向量 $v$ 选取 $v^*\dot v=0$ 的相位，对本征方程求导并在正交补上求逆，在谱隙大于一时有 $\|\dot v\|\le\|(I-vv^*)\dot H v\|$。同一荷的三个导数仍在正交的特征空间。先取其最大列范数，再对五荷用 Cauchy–Schwarz 与 $\sum q^2=66$，最后限制到所选的行，得到
$$
|s_S(x,\delta)-s_S(y,\delta')|\le\sqrt{132}\,d_1(x,y)+2\sqrt5\,|\delta-\delta'|.
$$
这里 $d_1$ 是两个短弧位移绝对值之和。在每个布局的一个平坦盲点 $z$ 周围，$d_1(x,z)<\varepsilon/\sqrt{132}$ 的菱形落在近盲区内。对 $0<\varepsilon\le1$，该菱形的归一化面积为 $\varepsilon^2/(264\pi^2)$。固定门槛可先把菱形略缩小，再用等分布、软化衰减和密度一满秩，得到同一下界。

对收缩门槛令 $L=\log T$，取 $J=\lfloor cL/\log L\rfloor$，使 $\|m\|_\infty\le J$ 时 $|u|\le L/(50\log L)$。三个相位及其逆乘 $1734$ 后均为代数整数。对 $e^{iu}-1$ 清分母取至多八次域范数，其余共轭的相位模均为一，得 $|u|\ge|e^{iu}-1|\ge2^{-7}1734^{-32\|m\|_\infty}$。引理 1 中最近素数幂的振荡积分绝对值至多 $T$，共振时等于 $T$。取 $x=e^{|u|}$，系数 $\Lambda(n_x)/\sqrt x$ 一致有界，其余误差除以 $\sqrt x$ 为 $T^{o(1)}$。引理 2 转成只含纵坐标的相位后，二者在截止范围内给出归一化模和界 $O(L^{-1}+\|m\|_\infty^2L^{-2})$。对 $\varepsilon$ 尺度下的固定光滑紧支撑测试函数，傅里叶系数绝对和、二阶矩、截止外尾和依次为 $O(1),O(\varepsilon^{-2}),O((\varepsilon J)^{-4})$，由二维格点逐层求和可得。采样与 Haar 均值之差是 $O(L^{-1}+\varepsilon^{-2}L^{-2}+(\varepsilon J)^{-4})=o(\varepsilon^2)$。[傅里叶转移的详细推导](https://github.com/the-omega-institute/trureturing/issues/35082#issuecomment-7232072134)

为在收缩菱形内也保留所有布局的严格满秩，取 $F$ 在 $z$ 处首个非零齐次 Taylor 项 $p_k$。在单位缩放后的菱形内去掉 $p_k=0$，损失面积为零。用内部紧集及支撑在其附近的光滑下逼近，可取得任意接近全菱形的面积，同时留有正的角度余量和 $|p_k|$ 下界。在这样的支撑上 $|F(z+\varepsilon y)|\ge c'\varepsilon^k$。将 $F$ 的各因子换为实际软化 Gram 行列式，紧致性与解析性给出乘积变化 $O(\delta_\gamma)$。在 $\sqrt T\le\gamma\le T$ 上它是 $O(T^{-1})=o(\varepsilon^k)$，较早零点仅占 $O(T^{-1/2})=o(\varepsilon^2)$。所选布局的奇异值由角度余量和软化界控制在 $\varepsilon$ 内，整个乘积仍非零。对光滑下逼近应用采样估计，再令其积分趋近菱形面积，得到所述下极限。

复算使用 Python 3.12、python-flint 0.9.0 及输入提交 `cb58c8e64ecea67184e0c1d68a690be098a7a384` 的 replay 目录。

```sh
python3 ch39-five-charge-zero-density.py /path/to/pinned/experiments/ch39-omega-research/replay
```

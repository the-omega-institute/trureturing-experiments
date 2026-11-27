# 凸片上外尔测度流对称部分模相位Hessian由曲率导数决定且纯超荷凸片上的射线定相参考部分以率零点一一三指数局域

测度流二阶变分的对称部分 $S$ 随参考测度的相位而变，变动为相位泛函的 Hessian；其中不随相位改变的部分是它的旋度，旋度恒等于曲率的一阶导数（在凸集或单连通区域上，对称部分模相位 Hessian 的等价类由旋度完全决定；全局情形未判定）。曲率导数有与体积无关的显式三链界：参考端点上，三条链的十五分量合计不超过 $7.929\times10^7\,e^{-0.07534622\,P}$，$P$ 为三个基点的周长，对一切相位选取成立。

在纯超荷切片上平凡场的凸邻域 $K_a=\{\max_g|A_g|\le1/4800\}$ 内，要求参考连接沿每条射线 $tA$ 的径向分量为零，就把参考测度的相位在 $K_a$（纯超荷切片）上的取值确定到只差一个常数。在这一定相规则下，对 $e,f$ 上的超荷坐标方向，参考连接只由终点投影决定，其对称部分不超过 $2.434\times10^{12}\,e^{-0.11301933\,d}$，$d\ge253$ 时小于一；加上 Kato 输运项，完整测度流的对称部分在任一正零点常数壳上不超过 $2.004\times10^{20}\,e^{-0.003047586\,d}$。弱方向与混合方向上的对称部分不受射线定相约束，仍可任意改变。

射线定相的规范不变性有一条判据：它在 $K_a$ 上规范不变，当且仅当十五分量加权的手征密度 $\mathcal A$ 在 $K_a$ 上恒为零。数值上（仅 $3^4$ 格），$\mathcal A$ 在场幅的四阶非零，指向射线定相不在容许参考类内，未证明。

本文没有给出的：这些界不是任何容许参考测度（包括 Kadoh–Kikukawa 与 Lüscher 的构造）的对称部分常数，“对所选构造的显式常数”仍未给出；射线定相能否延拓到 $K_a$ 以外，未判定。率 $0.11301933$ 只属于纯超荷凸片上的射线定相，不是模型实际所用参考测度的衰减率。

## 模型前提

沿用[外尔测度流的曲率部分显式指数局域而对称部分可随参考选择任意放大且衰减任意变慢](外尔测度流的曲率部分显式指数局域而对称部分可随参考选择任意放大且衰减任意变慢.md)与[十九点壳重叠核与测度输运项的显式指数局域性](十九点壳重叠核与测度输运项的显式指数局域性.md)的全部前提与记号：

- 周期四维格、周期 $\ell^1$ 距离 $d$、方格条件 $\|1-U_p\|\le\varepsilon=1/200$、一代十五个外尔分量与 $\sigma_c=\pm1$；
- 参考连接 $\mathcal L^{\rm ref}$、Kato 输运项 $j$、完整连接 $\mathcal L=\mathcal L^{\rm ref}+j$；
- 参考端点与各壳情形的 $\alpha,\beta,b,\rho_b,\ell$。

变分范数取十五个分量表示中像范数的最大值。

- **对称部分与曲率.** 对可交换的变分 $\delta_\eta,\delta_\zeta$，记 $S(\eta,\zeta)=\delta_\zeta\mathcal L_\eta+\delta_\eta\mathcal L_\zeta$，$F(\eta,\zeta)=\delta_\eta\mathcal L_\zeta-\delta_\zeta\mathcal L_\eta$。参考连接的 $F^{\rm ref}=\sum_c\sigma_cF_{P_0^c}$，完整连接的 $F=\sum_c\sigma_cF_{P_1^c}$。
- **可交换性.** 不同链上的左不变导数可交换。纯超荷切片上取坐标 $U_g=e^{iA_g}$（表示 $q$ 中为 $e^{iqA_g}$），坐标导数 $\partial_g$ 两两可交换，同一条链也可。$\partial_g$ 的像范数为 $\max|q|=6$。
- **凸片 $K_a$.** 纯超荷切片指弱链取恒等的场；$K_a=\{A:\max_g|A_g|\le a\}$，$a=1/4800$。
  - 在超荷 $q$ 表示中 $|\theta_p|\le4|q|a\le1/200$，故 $K_a$ 及其中每条射线 $\{tA\}$ 上的场都满足方格条件。
  - $A\mapsto e^{iA}$ 在 $K_a$ 上单射。$K_a$ 凸，含平凡场，在格点平移与超立方转动下不变，但不是规范不变的集合。
- **参考测度的定义域.** 参考测度只在文献的定义域内存在（体积足够大）。下文关于参考连接的结论都以参考测度在所论格上存在为前提；三链界本身只是关于投影的迹式，对任何格成立。

**新的模型选择（纯超荷凸片上的射线平行输运定相）.** 在 $K_a$ 上，把参考测度的相位取为使参考连接满足

$$
\mathcal L^{\rm rad}_A(A)=0\qquad(A\in K_a),
$$

即参考测度的基沿每条从平凡场出发的射线 $tA$ 平行输运。具体地，$\mathcal L^{\rm rad}=\mathcal L^{\rm ref}-i\,d\Phi$，$\Phi(A)=-i\int_0^1\mathcal L^{\rm ref}_{tA}(A)\,dt$（实值）。

- 这条规则只确定 $\Phi$ 在纯超荷切片 $K_a$ 上的取值；$\Phi$ 沿弱方向与混合方向的导数不受约束。
- 它不是 Kadoh–Kikukawa 或 Lüscher 构造本身的一部分。

## 结论

**定理 1（旋度恒等式）.** 对任何连接 $\mathcal L$ 与三条两两可交换的变分 $\delta_\theta,\delta_\eta,\delta_\zeta$，

$$
\delta_\theta F(\eta,\zeta)=\delta_\eta S(\zeta,\theta)-\delta_\zeta S(\eta,\theta).
$$

于是，若在某开集上 $S(\zeta,\theta)\equiv S(\eta,\theta)\equiv0$，则那里 $\delta_\theta F(\eta,\zeta)\equiv0$。相位重定义 $\mathcal L\mapsto\mathcal L-i\,d\Phi$ 不改变 $F$，因而不改变 $S$ 的旋度。

**定理 2（凸片上的完全分类）.** 在凸集 $K_a$ 上，下标与 $\eta,\zeta$ 都取 $K_a$ 的超荷坐标方向，$F$ 为参考连接的曲率。

1. 一个光滑的纯虚对称张量场 $T$ 是某个曲率为 $F$ 的纯虚连接（限制在超荷切片上）的对称部分，当且仅当 $\partial_jT_{ki}-\partial_kT_{ji}=\partial_iF_{jk}$ 对全部超荷坐标方向成立。
2. 两个曲率同为 $F$ 的这种连接相差 $-i\,d\Phi$，$\Phi$ 实值光滑；其对称部分相差 $-2i\,\partial_i\partial_j\Phi$。
3. 满足 $\mathcal L_A(A)=0$（$A\in K_a$）的这种连接恰有一个，即

$$
\mathcal L^{\rm rad}_A(\eta)=\int_0^1t\,F_{tA}(A,\eta)\,dt ,
$$

与之对应的相位 $\Phi$ 在 $K_a$ 上唯一到只差一个常数。因此在射线定相下，参考连接在超荷方向上只由终点三个标准重叠块的投影决定，与 Kadoh–Kikukawa 或 Lüscher 构造中相位的选取无关。它在格点平移与超立方转动下协变；反射不作断言。

$\mathcal L^{\rm rad}$ 的对称部分为

$$
S^{\rm rad}_A(\eta,\zeta)=\int_0^1t^2\bigl[(\partial_\zeta F)_{tA}(A,\eta)+(\partial_\eta F)_{tA}(A,\zeta)\bigr]\,dt ,
\qquad (\partial_\zeta F)(A,\eta)=\sum_gA_g\,\partial_\zeta F(\partial_g,\eta).
$$

弱方向与混合方向的对称部分不受射线定相约束，仍可按上述已有结果的定理 3 任意改变。

**定理 3（曲率导数的三链界）.** 设 $H$ 的谱含于 $[-\beta,-\alpha]\cup[\alpha,\beta]$，$b,\rho_b,\ell$ 取法同上，且 $\rho_b\ge1$。三条链 $g,e,f$ 的基点为 $z,x,y$，周长 $P=d(z,x)+d(x,y)+d(y,z)$。设三者两两可交换，且不全在同一条链上。则对十五个分量合计

$$
\Bigl|\sum_c\sigma_c\,\partial_\theta F_{P^c}(\eta,\zeta)\Bigr|\le C_3\,\|\theta\|\|\eta\|\|\zeta\|\,e^{-bP},\qquad
C_3=88\cdot k\cdot15\cdot\ell^3\rho_b^6e^{6b}.
$$

参考端点取 $k=2\cdot3$，十九点壳取 $k=2\cdot19$：

| 曲率导数 | 率 $b$（按周长） | $C_3$ | 界小于一 |
|---|---|---|---|
| 参考端点（三种壳情形共用） | $0.07534622$ | $(7.928,7.929)\times10^7$ | $P\ge242$ |
| 完整测度，任一正零点常数壳 | $0.001523793$ | $(1.512,1.513)\times10^{19}$ | $P\ge28983$ |

表中常数用有理上界 $e^{6b}\le1/(1-6b)$ 计算，并非最优；例如参考端点若用精确的 $e^{6b}$，$C_3$ 落在 $(6.82,6.84)\times10^7$ 内。由 $P\ge2\max\{d(x,y),d(y,z),d(z,x)\}$，任两条链相距 $d$ 时界不超过 $C_3e^{-2bd}$。由定理 1，这也是对称部分旋度的界，且对一切相位选取成立。

**定理 4（射线定相下对称部分的显式界）.** 在 $K_a$ 上、射线定相下，取两条不同的链 $e,f$（基点距离 $d$），$\eta,\zeta$ 为 $e,f$ 上的超荷坐标方向。对任意 $0<\lambda<1$，

$$
|S^{\rm rad}(\eta,\zeta)|\le\tfrac23\,a\cdot6\cdot C_3\cdot4\Bigl(\coth\tfrac{(1-\lambda)b}{2}\Bigr)^4\|\eta\|\|\zeta\|\,e^{-(1+\lambda)b\,d}.
$$

取参考端点的 $C_3,b$ 与 $\lambda=1/2$：

$$
|S^{\rm rad}(\eta,\zeta)|\le C_S\,\|\eta\|\|\zeta\|\,e^{-0.11301933\,d},\qquad C_S\in(2.433,2.434)\times10^{12},
$$

$d\ge253$ 时小于一。$C_S$ 用有理上界 $\coth(c/2)\le2(1+c)/c$ 与上述 $e^{6b}$ 上界计算，并非最优；仅换用精确的 $\coth$ 时落在 $(2.09,2.11)\times10^{12}$ 内。弱方向与混合方向的对称部分不受射线定相约束，仍可任意改变。

**推论（完整测度流对称部分）.** 在 $K_a$ 上、参考相位按射线定相，$\eta,\zeta$ 为 $e,f$ 上的超荷坐标方向。完整测度流的对称部分为 $S^{\rm rad}+S^j$，其中 $S^j(\eta,\zeta)=\delta_\zeta j_e(\eta)+\delta_\eta j_f(\zeta)$ 两项各由输运项定理 B 界定。于是

| 壳情形 | 率 | 前因子 | 其中输运项 | 界小于一 |
|---|---|---|---|---|
| 任一正零点常数壳 | $0.003047586$ | $(2.003,2.004)\times10^{20}$ | $2\times$ 定理 B | $d\ge15339$ |
| 11431 或 133931 常数壳 | $0.004442886$ | $<9.410\times10^{18}$ | $2\times$ 定理 B | $d\ge9834$ |
| 两者逐格任意切换 | $0.001700824$ | $<1.070\times10^{19}$ | $2\times$ 定理 B（两段之和） | $d\ge25762$ |

参考部分在前因子中所占比例小于 $10^{-7}$，且衰减快得多；这一界由 Kato 输运项主导。门槛 $n$ 按脚本中的精确前因子计：精确界在 $d=n$ 时小于一，在 $d=n-1$ 时不小于一。弱方向与混合方向的对称部分不受射线定相约束，仍可任意改变。

**命题 5（射线定相的规范判据）.** 对格点实函数 $\omega$，记规范方向 $\nabla\omega$ 为超荷坐标中的常向量场 $(\nabla\omega)_{(x,\mu)}=\omega(x+\hat\mu)-\omega(x)$。令

$$
\mathcal A(x)(A)=\sum_c\sigma_cq_c\,\operatorname{tr}P_{0,+}^{(c)}(A)(x,x)-\bigl(\text{同式在 }A=0\text{ 的值}\bigr),
$$

迹取 $x$ 处的自旋、内部与表示指标，$P_{0,+}^{(c)}$ 为分量 $c$ 的终点正谱投影，$q_c$ 为其超荷。在 $K_a$ 的内部：

1. $F^{\rm ref}_A(\nabla\omega,\eta)=i\sum_x\omega(x)\,\delta_\eta\mathcal A(x)(A)$；
2. $\mathcal L^{\rm rad}_A(\nabla\omega)=-i\sum_x\omega(x)\bigl[\mathcal A(x)(A)-\bar{\mathcal A}(x)(A)\bigr]$，其中 $\bar{\mathcal A}(A)=\int_0^1\mathcal A(tA)\,dt$；
3. $\delta_{\nabla\omega}\mathcal L^{\rm rad}_\eta=i\sum_x\omega(x)\,\delta_\eta\bar{\mathcal A}(x)$。

因此，射线定相的测度流沿全部规范方向不变（对一切 $\omega$ 与超荷方向 $\eta$，$\delta_{\nabla\omega}\mathcal L^{\rm rad}_\eta=0$），当且仅当 $\mathcal A\equiv0$ 于 $K_a$。

**命题 6（由射线外的一点判定）.** 设 $A$ 为纯超荷场（不必在 $K_a$ 内），各分量的终点核 $H(tA)$ 对全部 $t\in[0,1]$ 有隙。若某个 $x$ 处 $\mathcal A(x)(A)\ne0$，则 $\mathcal A$ 在 $K_a$ 上不恒为零，从而射线定相在 $K_a$ 上不规范不变，不在容许参考类内。

**数值（仅 $3^4$ 格）.** 取 Wilson 质量 $-1$ 的标准重叠核、十五分量权与一个固定随机场 $T$，场 $\varepsilon T$：

- $\max_x|\mathcal A|$ 在 $\varepsilon=10^{-2}$、$5\times10^{-3}$、$2.5\times10^{-3}$ 时分别落在 $(2.10,2.11)\times10^{-6}$、$(1.32,1.33)\times10^{-7}$、$(8.26,8.28)\times10^{-9}$ 内；
- 相邻两级之比都在 $(15.5,16.5)$ 内，即 $\varepsilon^4$ 行为；
- $|\sum_x\mathcal A|$ 小于 $10^{-12}$；随机规范变换前后之差小于 $10^{-12}$；零场值小于 $10^{-14}$；
- 射线 $\{t\cdot10^{-2}T\}$ 上的浮点隙下界大于 $0.9$：在 $21$ 个等距格点上取浮点最小 $|\lambda|$，相邻格点间用 Weyl 不等式与 $\|H(t)-H(s)\|\le8|q|\max|A|\,|t-s|$ 插值；
- 最大值至少是 $|\sum_x\mathcal A|$、规范变换差、零场值与 $10^{-12}$ 中最大者的一千倍。

这些场幅的 $\max|A_g|$ 都大于 $a=1/4800$，即场不在 $K_a$ 内；由命题 6，只要射线上有隙，$K_a$ 外一点的非零值就足以判定。由于本征分解是浮点计算，结论只是数值：$\mathcal A$ 在四阶非零，指向射线定相不在容许参考类内，未证明。$3^4$ 是否属于参考测度的体积定义域，也未判定。

## 证明

**定理 1.** 由 $S(\zeta,\theta)=\delta_\theta\mathcal L_\zeta+\delta_\zeta\mathcal L_\theta$ 与混合导数可交换，

$$
\delta_\eta S(\zeta,\theta)-\delta_\zeta S(\eta,\theta)=\delta_\eta\delta_\theta\mathcal L_\zeta-\delta_\zeta\delta_\theta\mathcal L_\eta=\delta_\theta\bigl(\delta_\eta\mathcal L_\zeta-\delta_\zeta\mathcal L_\eta\bigr)=\delta_\theta F(\eta,\zeta).
$$

**定理 2.**

- **第 2 条.** 两连接之差是 $K_a$ 上的闭纯虚一形式。由 Poincaré 引理，它等于 $-i\,d\Phi$，$\Phi$ 实值；对称部分之差为 $-2i\,\partial_i\partial_j\Phi$。这一步用到 $K_a$ 凸。在不单连通的区域上，闭形式未必恰当，所以这里的结论不推广到全局场空间。
- **第 1 条.** 必要性即定理 1。充分性：任取曲率为 $F$ 的连接 $\mathcal L^0$（例如 $\mathcal L^{\rm ref}$），其对称部分 $S^0$ 满足同一旋度方程，故 $D=T-S^0$ 满足 $\partial_jD_{ki}=\partial_kD_{ji}$。对固定的 $i$，$D_{\cdot i}$ 是闭一形式，$D_{ki}=\partial_k\psi_i$。又 $D$ 对称，$\partial_k\psi_i=\partial_i\psi_k$，故 $\psi=d\chi$，$D=\partial\partial\chi$，$\chi$ 纯虚。于是 $\mathcal L^0+\tfrac12d\chi$ 的曲率仍为 $F$，对称部分为 $T$。
- **第 3 条.** 令 $\varphi(A)=\int_0^1\mathcal L^{\rm ref}_{tA}(A)\,dt$。由 $\frac{d}{dt}\bigl[t\,\mathcal L^{\rm ref}_{tA}(\eta)\bigr]=\mathcal L^{\rm ref}_{tA}(\eta)+t\,(\partial_A\mathcal L^{\rm ref}_\eta)_{tA}$ 在 $[0,1]$ 上积分，

$$
\mathcal L^{\rm ref}_A(\eta)-\partial_\eta\varphi(A)=\int_0^1t\bigl[(\partial_A\mathcal L^{\rm ref}_\eta)-(\partial_\eta\mathcal L^{\rm ref}_A)\bigr]_{tA}dt=\int_0^1t\,F_{tA}(A,\eta)\,dt .
$$

  这就是 $\mathcal L^{\rm rad}=\mathcal L^{\rm ref}-d\varphi$，$\Phi=-i\varphi$。取 $\eta=A$ 得 $\mathcal L^{\rm rad}_A(A)=0$。若另有 $\mathcal L^{\rm rad}-i\,d\Psi$ 也满足径向条件，则 $A\cdot\nabla\Psi(A)=0$，$\Psi$ 沿每条射线为常数；由在 $0$ 处连续，$\Psi$ 为常数。
- **协变性.** 平移与超立方转动在 $A$ 坐标上线性作用，把 $K_a$ 映为自身、固定平凡场并保持射线。重叠投影在其下协变：转动的自旋矩阵与 $\gamma_5$ 可交换。所以 $F^{\rm ref}$ 协变，$\mathcal L^{\rm rad}$ 由 $F^{\rm ref}$ 经同一公式给出，也协变。反射把 $\gamma_5$ 变号，交换手征，此处不作断言。
- **对称部分.** $\partial_\zeta\mathcal L^{\rm rad}_A(\eta)=\int_0^1t\bigl[F_{tA}(\zeta,\eta)+t\,(\partial_\zeta F)_{tA}(A,\eta)\bigr]dt$。对 $\eta,\zeta$ 对称化后 $F(\zeta,\eta)+F(\eta,\zeta)=0$，余下即所写公式。
- **适用方向.** 全部计算都在超荷坐标中进行，只用到 $\Phi$ 在 $K_a$ 上的值。$\Phi$ 沿弱方向的导数不出现在径向条件中，所以弱方向与混合方向的对称部分不受约束。

**定理 3.** 写 $F_P(\eta,\zeta)=\operatorname{Tr}(P\,\delta_\eta P\,\delta_\zeta P)-\operatorname{Tr}(P\,\delta_\zeta P\,\delta_\eta P)$，$P=\oint R$，$\delta_\eta P=\oint RE_\eta R$，

$$
\delta_\theta\delta_\eta P=\oint\bigl(RE_\theta RE_\eta R+RE_\eta RE_\theta R+RE_{\theta\eta}R\bigr),\qquad E_{\theta\eta}=\delta_\theta\delta_\eta H .
$$

$H$ 对每条链的链变量分别线性，所以不同链上 $E_{\theta\eta}=0$。同一条链上（只在纯超荷坐标中出现），$E_{\theta\eta}$ 与 $E_\theta$ 一样只有两块，范数至多 $\|\theta\|\|\eta\|$（表示 $q$ 中为 $q^2$）。

- **展开.** $\delta_\theta\operatorname{Tr}(P\,\delta_\eta P\,\delta_\zeta P)$ 展开为三部分：$\operatorname{Tr}(\delta_\theta P\,\delta_\eta P\,\delta_\zeta P)$ 一项；$\operatorname{Tr}(P\,\delta_\theta\delta_\eta P\,\delta_\zeta P)$ 两项，加一项含 $E_{\theta\eta}$；$\operatorname{Tr}(P\,\delta_\eta P\,\delta_\theta\delta_\zeta P)$ 两项，加一项含 $E_{\theta\zeta}$。另一次序相同。
- **计数.** 合计十项，各含三个围道积分、六个 $R$ 与 $E_\theta,E_\eta,E_\zeta$ 各一次。三条链不全重合时，$E_{\theta\eta}$ 与 $E_{\theta\zeta}$ 至多一个非零，故另有至多两项，各含三个围道积分、五个 $R$ 与两个插入。若三条链全重合，纯超荷中 $\eta\propto\zeta$，$F(\eta,\zeta)\equiv0$。
- **六预解式项.** 可循环写成 $\operatorname{Tr}(E_\theta XE_\eta YE_\zeta Z)$，$N_b(X)N_b(Y)N_b(Z)\le\rho_b^6$。每个 $E$ 两块，展开为八组；每组块的秩至多 $kr$。三段核分别从一块的端点连到下一块的端点，每段距离至少为相应基点距离减二。所以每组不超过 $kr\,\|\theta\|\|\eta\|\|\zeta\|\rho_b^6e^{-b(P-6)}$。
- **五预解式项.** 含两个插入、四组块，周长 $P=2d$（两条链重合）。每组不超过 $kr\,\|\theta\|\|\eta\|\|\zeta\|\rho_b^5e^{-b(P-4)}\le kr\,\|\theta\|\|\eta\|\|\zeta\|\rho_b^6e^{-b(P-6)}$，用到 $\rho_b\ge1$。

于是常数为 $(10\cdot8+2\cdot4)\,kr\,\ell^3\rho_b^6e^{6b}$。对十五个分量求和时，以 $\sum_cr_c\|\theta_c\|\|\eta_c\|\|\zeta_c\|\le15\|\theta\|\|\eta\|\|\zeta\|$ 代替 $r\|\theta\|\|\eta\|\|\zeta\|$。参考端点上只有三个标准块的投影随场变化，取 $k=6$；谱参数与 $b,\rho_b,\ell$ 取自曲率局域性的证书，与零点无关。估计不涉及格点求和，常数与体积无关。

**定理 4.** 射线上每点都在 $K_a$ 内，满足方格条件，所以定理 3 的常数一致成立。由定理 2 的公式与 $|A_g|\le a$、$\|\partial_g\|=6$，

$$
|S^{\rm rad}(\eta,\zeta)|\le\int_0^1t^2dt\cdot2\cdot a\cdot6\cdot C_3\|\eta\|\|\zeta\|\sum_ge^{-bP(g)} .
$$

链 $g=(z,\mu)$ 每个格点四条，$P(g)=d(x,y)+d(z,x)+d(z,y)$。由 $d(z,x)+d(z,y)\ge d(x,y)$，

$$
\sum_ge^{-bP(g)}\le4e^{-(1+\lambda)b\,d(x,y)}\sum_ze^{-(1-\lambda)b\,d(z,x)} .
$$

把周期格上每点取各坐标在 $(-L/2,L/2]$ 内的代表，得到到 $\mathbb Z^4$ 的单射，且保持到 $x$ 的周期 $\ell^1$ 距离。故末一和不超过 $\sum_{n\in\mathbb Z^4}e^{-c|n|_1}=\coth(c/2)^4$，$c=(1-\lambda)b$。脚本用 $e^{-c}\le1/(1+c)$ 得有理上界 $\coth(c/2)\le2(1+c)/c$。

**推论.** $\mathcal L=\mathcal L^{\rm rad}+j$ 的对称部分为 $S^{\rm rad}+S^j$。输运项定理 B 对任一对不同链给出 $|\delta_\zeta j_e(\eta)|\le C_B\|\eta\|\|\zeta\|e^{-2bd}$，两项相加得 $2C_B$；$K_a$ 满足其方格条件。参考部分的率 $0.11301933$ 大于各壳情形的率，可并入较慢的率。

**命题 5.**

- **规范方向上的投影变分.** 沿 $+\nabla\omega$ 移动场，等于施加规范变换 $g=e^{-iq\omega}$：链 $U(x,\mu)\mapsto g(x)U(x,\mu)g(x+\hat\mu)^{-1}$ 使 $A_{(x,\mu)}$ 增加 $\omega(x+\hat\mu)-\omega(x)$。投影协变，$P(gU)=gP(U)g^{-1}$，故 $\delta_{\nabla\omega}P=[X,P]$，$X=-iq\omega$（逐站乘子）。
- **一条恒等式.** 记 $Q=\delta_\eta P$，它满足 $PQP=0$ 与 $Q=PQ+QP$。展开

$$
\operatorname{Tr}(P[X,P]Q)-\operatorname{Tr}(PQ[X,P])=\operatorname{Tr}(PXPQ)-\operatorname{Tr}(PXQ)-\operatorname{Tr}(PQXP)+\operatorname{Tr}(PQPX).
$$

  其中第一项与第四项为零，第二、三项之和为 $\operatorname{Tr}(X(QP+PQ))$。所以 $F_P(\nabla\omega,\eta)=-\operatorname{Tr}(X\,\delta_\eta P)=iq\sum_x\omega(x)\,\delta_\eta\operatorname{tr}P(x,x)$。
- **第 1 条.** 对十五个分量乘 $\sigma_c$ 求和：右手分量取 $1-P_+$，按约定 $F_{1-P_+}=-F_{P_+}$ 已计入 $\sigma_c$。十六个常投影块的迹与场无关，被 $\delta_\eta$ 消去；减去零场值也不影响 $\delta_\eta$。这就得到第 1 条。
- **第 2 条.** 由定理 2 第 3 条，

$$
\mathcal L^{\rm rad}_A(\nabla\omega)=-\int_0^1t\,F_{tA}(\nabla\omega,A)\,dt=-i\sum_x\omega(x)\int_0^1t\,\frac{d}{dt}\mathcal A(x)(tA)\,dt .
$$

  分部积分，$\int_0^1t\,\tfrac{d}{dt}\mathcal A(tA)\,dt=\mathcal A(A)-\bar{\mathcal A}(A)$，即第 2 条。
- **第 3 条.** $\nabla\omega$ 与超荷坐标方向都是常向量场，可交换，所以 $\delta_{\nabla\omega}\mathcal L_\eta=F(\nabla\omega,\eta)+\delta_\eta\mathcal L_{\nabla\omega}$。代入第 1、2 条，$\delta_\eta\mathcal A$ 两项相消，得第 3 条。
- **判据.** 若 $\mathcal A\equiv0$ 于 $K_a$，则 $\bar{\mathcal A}\equiv0$，规范不变成立。反之，若对一切 $\omega,\eta$ 有 $\delta_{\nabla\omega}\mathcal L^{\rm rad}_\eta=0$，则取 $\omega$ 为单点的指示函数，得 $\delta_\eta\bar{\mathcal A}(x)=0$ 对一切 $x$ 与超荷方向 $\eta$ 成立。于是 $\bar{\mathcal A}(x)$ 在连通的 $K_a$ 内部为常数，等于零场值 $0$。对 $s\in(0,1]$，$s\bar{\mathcal A}(sA)=\int_0^s\mathcal A(uA)\,du=0$，对 $s$ 求导得 $\mathcal A(sA)=0$。由连续性，$\mathcal A\equiv0$ 于 $K_a$。

**命题 6.**

- **解析性.** $H(tA)$ 对 $t$ 是整函数（链变量为 $e^{iqtA_g}$）。$[0,1]$ 紧且每点有隙，所以存在复邻域，在其上 Riesz 投影 $\oint(z-H(tA))^{-1}dz$ 解析。因而 $t\mapsto\mathcal A(x)(tA)$ 在含 $[0,1]$ 的开区间上实解析。
- **反证.** 若 $\mathcal A\equiv0$ 于 $K_a$，则对 $0\le t\le a/\max_g|A_g|$ 有 $\mathcal A(tA)=0$。由解析延拓，$\mathcal A(tA)$ 在 $[0,1]$ 上恒为零，与 $\mathcal A(x)(A)\ne0$ 矛盾。
- **结论.** 由命题 5，射线定相不规范不变，因而不满足容许参考类的规范不变性一条。

## 关键数字

- **凸片.** $a=1/4800$，$4\cdot6\cdot a=1/200$。
- **三链前因子.** 参考端点 $C_3\in(7.928,7.929)\times10^7$；普适壳 $C_3\in(1.512,1.513)\times10^{19}$。$88=10\cdot8+2\cdot4$。均为有理上界，非最优。
- **射线定相.** 率 $\tfrac32b_0=0.11301933$，$C_S\in(2.433,2.434)\times10^{12}$，门槛 $253$。有理上界，非最优。
- **完整测度对称部分.** 三种壳情形的前因子与门槛见推论表；普适情形前因子在 $(2.003,2.004)\times10^{20}$ 内。参考部分占比小于 $10^{-7}$。
- **门槛.** 对精确有理前因子 $C$ 与率 $r$，以 Arb 核对 $Ce^{-rn}<1\le Ce^{-r(n-1)}$。
- **手征密度（仅数值，$3^4$）.** 见命题 6 之后的数值条目。

## 与已有结果的关系

[外尔测度流的曲率部分显式指数局域而对称部分可随参考选择任意放大且衰减任意变慢](外尔测度流的曲率部分显式指数局域而对称部分可随参考选择任意放大且衰减任意变慢.md)的定理 3 证明对称部分不由投影与容许性决定，并写道“只能对所选的具体参考构造另行估计；本文未对这些具体构造给出常数”。本文推进这一句的三个方面：

- **相位无关的内容.** 定理 1、2 说明，在凸片上，对称部分中与相位无关的内容恰是旋度，即曲率导数；其余部分恰为相位泛函的 Hessian，正是该文定理 3 所显式化的任意性。
- **显式常数.** 定理 3 给出旋度的显式常数，对一切相位选取成立。
- **定相规则.** 定理 4 给出一条在纯超荷凸片上把相位确定到常数的定相规则及其显式界。命题 5、6 给出它的规范不变性判据，数值指向不成立。

那句话本身仍成立：本文没有估计 Kadoh–Kikukawa 或 Lüscher 原有相位下的对称部分常数，也没有给出任何容许参考测度的对称部分常数。

该文推论说参考连接在曲率非零处不能经相位重定义消去。定理 1 给出对称部分的同类结论：凡 $\partial_\theta F^{\rm ref}(\eta,\zeta)\not\equiv0$ 的开集上，$S(\zeta,\theta)$ 与 $S(\eta,\theta)$ 不能同时恒为零。

推论中的输运项部分直接取自[十九点壳重叠核与测度输运项的显式指数局域性](十九点壳重叠核与测度输运项的显式指数局域性.md)定理 B 的精确前因子，本文只把两项相加。

## 证据层级

- **书面证明.** 定理 1–4、推论、命题 5 与 6。
  - 包括十项加两项的计数、八组与四组块、$e^{6b}$ 因子，以及格点求和的单射论证。
  - 包括命题 5 的恒等式 $F_P(\nabla\omega,\eta)=-\operatorname{Tr}(X\delta_\eta P)$ 与当且仅当判据，以及命题 6 的解析延拓。
  - 脚本的 `three_link_constant`、`ray_symmetric_constant` 按这些公式计算常数，不复算推导本身。
- **精确有理证书.**
  - $K_a$ 的方格条件；
  - 三链前因子、射线定相前因子与率；
  - 完整测度对称部分的三个前因子与率。

  谱参数取自重跑的[曲率局域性脚本](ch39-shell-measure-curvature-locality.py)及其调用的输运项证书。
- **Arb 区间证书.** 全部距离与周长门槛。
- **仅数值.**
  - *玩具族.* 在 8 维、秩 4 的三参数投影族上，由一个光滑基的连接作有限差分：
    - 旋度恒等式两端相对差小于 $2\times10^{-3}$；
    - 射线定相连接的径向分量小于 $10^{-5}$；
    - $\mathcal L^{\rm rad}$ 与 $\int_0^1tF_{tA}(A,\eta)dt$ 之差小于 $10^{-5}$；
    - $S^{\rm rad}$ 与公式的相对差小于 $2\times10^{-3}$。
  - *$3^4$ 格曲率导数.* 取 Wilson 质量 $-1$ 的标准重叠核、一个满足 $|q|=6$ 方格条件的随机小 $U(1)$ 场（幅 $10^{-4}$，该场未必在 $K_a$ 内）与三条不同的链 $e,f,g$。
    - 十五分量端点曲率对 $g$ 的中心差分落在 $(2.2,2.3)\times10^{-9}$ 内（乘 $i$）。
    - 步长减半的相对差小于 $10^{-2}$；场幅减半时比值在 $(3.9,4.1)$ 内，与二次标度一致。
    - 零场处同一差分小于其百分之一。
  - *$3^4$ 格手征密度 $\mathcal A$.* 见命题 6 之后的数值条目。场不在 $K_a$ 内；隙与 $\mathcal A$ 都由浮点本征分解得到。
  - $3^4$ 是否属于参考测度的体积定义域，未判定。

## 复算

    python3 experiments/ch39-omega-research/replay/ch39-measure-symmetric-part-ray-gauge.py

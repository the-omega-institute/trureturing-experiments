# 纯超荷凸片内十五分量手征密度在二的四次方格上有Arb证书非零因而该格上射线定相不规范不变

在 $2^4$ 周期格上取纯超荷场 $A_*=T/5000$（$T_g=\pm1$，逐链给定，$\max_g|A_{*,g}|=1/5000<1/4800$，故 $A_*$ 在凸片 $K_a$ 的内部）。十五分量加权的手征密度在格点 $0$ 处落在 $(2.0362339,2.0362340)\times10^{-11}$ 内，在十六个格点上全部非零，绝对值都大于 $2.44\times10^{-12}$。在 $2\times2\times1\times1$ 格（四个格点）与 $3^4$ 格上，同样有 $K_a$ 内部的场使它非零：分别落在 $(-2.6783766,-2.6783765)\times10^{-11}$ 与 $(-9.0757235,-9.0757234)\times10^{-13}$ 内。

这里的射线定相直接定义为 $K_a$ 上由终点投影的曲率给出的连接 $\mathcal L^{\rm rad}_A(\eta)=\int_0^1tF_{tA}(A,\eta)\,dt$（引理 0），不经过参考测度。由[凸片上外尔测度流对称部分模相位Hessian由曲率导数决定且纯超荷凸片上的射线定相参考部分以率零点一一三指数局域](凸片上外尔测度流对称部分模相位Hessian由曲率导数决定且纯超荷凸片上的射线定相参考部分以率零点一一三指数局域.md)的命题 5，在这三个格上，这个连接沿规范方向不是不变的，所以它不满足容许参考类的规范不变性一条。$K_a$ 上处处有隙，由初等界 $|\operatorname{spec}H_q|\ge199/200$ 给出（引理 1）。手征密度的奇数阶项因电荷共轭对称全部为零，首项为 $12960\,f_4$（定理 3）。

在超立方格 $L^4$ 中，$L=2$ 是最小者：$L=1$ 时没有非平凡规范方向，且 $\mathcal A\equiv0$。在长方格中，$2\times2\times1\times1$（体积 $4$）有证书。体积 $2$、$3$ 的长方格（$2\times1\times1\times1$、$3\times1\times1\times1$ 及其坐标置换）上，各一个随机场的浮点读数 $\max_x|\mathcal A|/3$ 小于 $10^{-13}$（仅数值）；$\mathcal A$ 在那里是否恒为零，未判定。

本文没有给出的：这三个格是否属于参考测度的体积定义域，未判定。所引 Lüscher 定理 5.4 的陈述只要求 $L/\varrho\ge n$，未给出 $n$ 的数值，由此无法判定这三个格是否满足这一条件。因此，在参考测度定义域内的体积上射线定相是否同样不规范不变，仍未判定。

## 模型前提

沿用[凸片上外尔测度流对称部分模相位Hessian由曲率导数决定且纯超荷凸片上的射线定相参考部分以率零点一一三指数局域](凸片上外尔测度流对称部分模相位Hessian由曲率导数决定且纯超荷凸片上的射线定相参考部分以率零点一一三指数局域.md)（下称“射线定相一文”）与[外尔测度流的曲率部分显式指数局域而对称部分可随参考选择任意放大且衰减任意变慢](外尔测度流的曲率部分显式指数局域而对称部分可随参考选择任意放大且衰减任意变慢.md)的全部前提与记号。

- **格与场.** 周期长方格 $\Lambda=\mathbb Z_{L_0}\times\cdots\times\mathbb Z_{L_3}$。纯超荷场 $A_g\in\mathbb Q$，超荷 $q$ 表示中的链变量为 $e^{iqA_g}$。$K_a=\{\max_g|A_g|\le a\}$，$a=1/4800$。
- **终点核.** $H_q(A)=\gamma_5D_W$，$D_W$ 为 Wilson 质量 $-1$ 的标准 Wilson 算子：

$$
D_W=3-\sum_\mu\bigl(P^-_\mu\otimes T_\mu^{U}+P^+_\mu\otimes(T_\mu^{U})^\dagger\bigr),\qquad P^\pm_\mu=\tfrac12(1\pm\gamma_\mu),
$$

  $T_\mu^U$ 为带链变量的前向平移。$\gamma$ 矩阵取射线定相一文数值脚本的手征基：$\gamma_k=\begin{pmatrix}0&-i\sigma_k\\ i\sigma_k&0\end{pmatrix}$（$k=1,2,3$），$\gamma_4=\begin{pmatrix}0&1\\1&0\end{pmatrix}$，$\gamma_5=\gamma_1\gamma_2\gamma_3\gamma_4=\operatorname{diag}(1,1,-1,-1)$。$P_+=(1+\operatorname{sign}H)/2$。
- **十五分量手征密度.** 与射线定相一文命题 5 相同：

$$
\mathcal A(x)(A)=\sum_c\sigma_cq_c\operatorname{tr}P_{0,+}^{(c)}(A)(x,x)-\bigl(\text{同式在 }A=0\text{ 的值}\bigr).
$$

  - 迹取自旋、内部与表示指标。每个分量的终点在内部十九维上含三个相同的标准块与十六个常投影块，常投影块的迹被零场值消去。
  - 按超荷归并，权 $w_q=\sigma\cdot$重数为 $w_1=6$，$w_{-3}=2$，$w_4=-3$，$w_{-2}=-3$，$w_{-6}=-1$（十五个分量）。于是

$$
\mathcal A(x)(A)=3\sum_qw_q\,q\,\bigl[\operatorname{tr}_{\rm spin}P_+^{(q)}(A)(x,x)-\operatorname{tr}_{\rm spin}P_+^{(q)}(0)(x,x)\bigr].
$$

- **三个对象的区分.**
  - *模型实际所用的参考测度*：纯超荷扇区取 Lüscher 的构造，它在其定义域内按构造规范不变。
  - *射线定相*：本文按引理 0 定义的 $K_a$ 上的连接 $\mathcal L^{\rm rad}$。在参考测度存在的格上，它与射线定相一文的定义（参考测度的一种相位选取）一致。
  - *文献构造*：Lüscher（[hep-lat/9811032](https://arxiv.org/abs/hep-lat/9811032)）与 Kadoh–Kikukawa（[arXiv:0709.3658](https://arxiv.org/abs/0709.3658)）。

  本文的结论只关于射线定相。
- **容许参考类的规范不变性.** 与射线定相一文相同：对一切格点函数 $\omega$ 与超荷坐标方向 $\eta$，$\delta_{\nabla\omega}\mathcal L_\eta=0$。

## 结论

**引理 0（不经参考测度的射线定相）.** 设 $K_a$ 上各 $H_q$ 有隙（引理 1）。令 $F=\sum_c\sigma_cF_{P_0^c}$，$F_P=\operatorname{Tr}(P\,dP\wedge dP)$，即 $F_P(\eta,\zeta)=\operatorname{Tr}(P[\delta_\eta P,\delta_\zeta P])$，它只由终点投影决定。在超荷坐标中定义

$$
\mathcal L^{\rm rad}_A(\eta)=\int_0^1t\,F_{tA}(A,\eta)\,dt\qquad(A\in K_a).
$$

则 $F$ 是 $K_a$ 上的闭二形式，$\mathcal L^{\rm rad}$ 的曲率 $d\mathcal L^{\rm rad}$ 等于 $F$，且 $\mathcal L^{\rm rad}_A(A)=0$。射线定相一文命题 5 第 1–3 条及其判据的证明只用到这三点、$F$ 的规范协变性与 $K_a$ 的凸性，因而对这里定义的 $\mathcal L^{\rm rad}$ 在任何格上成立，不需要参考测度存在。

**引理 1（$K_a$ 上的一致隙）.** 对任何周期长方格、任何 $A\in K_a$ 与 $|q|\le6$，

$$
|\operatorname{spec}H_q(A)|\ge1-4\cdot6\cdot a=\tfrac{199}{200}.
$$

**定理 2（区间证书）.** 记 $T^{(\Lambda)}$ 为脚本中逐字给出的 $\pm1$ 模式，$A_*=T^{(\Lambda)}/5000$。

| 格 $\Lambda$ | 格点 $x$ | $\mathcal A(x)(A_*)$ | 最大方格角 | 五种超荷上 $\min\lvert\operatorname{spec}H_q(A_*)\rvert$ 的下界 |
|---|---|---|---|---|
| $2^4$ | $0$ | $(2.0362339,2.0362340)\times10^{-11}$ | $1/2500$ | $0.999100$ |
| $2\times2\times1\times1$ | $0$ | $(-2.6783766,-2.6783765)\times10^{-11}$ | $1/2500$ | $0.998330$ |
| $3^4$ | $30$，即 $(1,0,1,0)$ | $(-9.0757235,-9.0757234)\times10^{-13}$ | $1/1250$ | $0.999086$ |

- 三个场都满足 $\max_g|A_{*,g}|=1/5000<a$，都在 $K_a$ 的内部。超荷 $6$ 下的方格角不超过 $6/1250<1/200$。
- 表中隙下界由 Bauer–Fike 界认证，都不小于引理 1 的 $199/200$，与之相容。
- 在 $2^4$ 上，十六个格点处 $\mathcal A(x)(A_*)$ 都不含零，绝对值都大于 $2.44\times10^{-12}$。

**推论（这三个格上射线定相不规范不变）.** 在 $2^4$、$2\times2\times1\times1$、$3^4$ 上，$\mathcal A$ 在 $K_a$ 上不恒为零。由引理 0 与射线定相一文命题 5，存在 $\omega$ 与超荷方向 $\eta$，使 $\delta_{\nabla\omega}\mathcal L^{\rm rad}_\eta\ne0$ 在 $K_a$ 内部某处成立。因此在这三个格上，射线定相不满足容许参考类的规范不变性。由于 $A_*\in K_a$，不需要该文命题 6 的解析延拓。

**命题 2（总和恒为零）.** 在 $K_a$ 上，对任何格，$\sum_x\mathcal A(x)(A)\equiv0$。于是在 $1^4$ 上 $\mathcal A\equiv0$。

**定理 3（奇数阶为零与首项）.** 记 $f(A)(x)=\operatorname{tr}_{\rm spin}P_+^{(1)}(A)(x,x)$。

1. $H_q(A)=H_1(qA)$。在 $0$ 的某邻域内，$f$ 是 $A$ 的收敛幂级数 $f=\sum_kf_k$，$f_k$ 为 $k$ 次齐次多项式。作为幂级数系数的恒等式，$f(qA)=\sum_kq^kf_k(A)$，从而

$$
\mathcal A(x)(A)=3\sum_{k\ge1}M_{k+1}\,f_k(A)(x),\qquad M_m=\sum_qw_qq^m ,
$$

   其中 $M_1=0$，$M_2=-72$，$M_3=0$，$M_4=-1944$，$M_5=4320$。
2. （电荷共轭）在 $K_a$ 上 $f(-A)=f(A)$，故一切奇数阶 $f_k$ 为零。
3. 于是 $\mathcal A=3\sum_{k\ \rm even}M_{k+1}f_k$。二阶项因 $M_3=0$ 为零，所以

$$
\mathcal A(x)(A)=12960\,f_4(A)(x)+O(|A|^6).
$$

**命题 4（最小体积）.**

- 在超立方格 $L^4$ 中，射线定相在 $K_a$ 上不规范不变的最小 $L$ 是 $2$。
  - $L=1$ 时 $\nabla\omega=0$，规范方向只有零向量，规范不变性自动成立；此外由命题 2，$\mathcal A\equiv0$。
  - $L=2$ 时由推论不成立。
- 在长方格中，体积 $4$ 的 $2\times2\times1\times1$ 已有证书，体积 $1$ 的情形同上。
- 体积 $2$、$3$ 的长方格只有 $2\times1\times1\times1$、$3\times1\times1\times1$ 及其坐标置换。由超立方转动（某一平面内的 $90^\circ$ 转动，配以与 $\gamma_5$ 可交换的自旋转动）的协变性，坐标置换后的格上 $\mathcal A$ 是否恒为零与原格相同。
  - 脚本只测了 $2\times1\times1\times1$、$3\times1\times1\times1$、$4\times1\times1\times1$ 三种形状。
  - 对每一种取一个随机高斯场 $T$，在 $\varepsilon T$（$\varepsilon=2\times10^{-2},10^{-2}$）上，浮点读数 $\max_x|\mathcal A|/3$ 小于 $10^{-13}$。
  - 这只是数值；$\mathcal A$ 在这些格上是否恒为零，未判定。
  - 所以在长方格中，最小体积是 $2$、$3$、$4$ 中的哪一个，未判定。

**定义域（未判定）.** Lüscher 构造在有限体积上的存在性来自其定理 5.4。按 arXiv 版 hep-lat/9811032 的 PDF 文本，定理 5.4 写道：“If the lattice is sufficiently large compared to the localization range $\varrho$ of the Dirac operator, say $L/\varrho\ge n$”，又写明 $n,\kappa,\nu$ 与场无关；定理 5.4 的陈述中未给出它们的数值。原文的容许条件为 $|F_{\mu\nu}(x)|<\epsilon$，$|e_\alpha|\epsilon<1/30$，本文的三个场满足它：$6\cdot\tfrac1{1250}<\tfrac1{30}$。由此无法判定 $2^4$、$2\times2\times1\times1$、$3^4$ 是否满足 $L/\varrho\ge n$。

## 证明

**引理 0.**

- **闭性.** 由 $P^2=P$ 得 $dP=P\,dP+dP\,P$，于是 $P\,dP\,P=0$ 且 $P^\perp dP\,P^\perp=0$，$P^\perp=1-P$。所以 $dP$ 把 $\operatorname{ran}P$ 映入 $\operatorname{ran}P^\perp$，反之亦然。三个这样的“交换两子空间”的算子之积仍交换两子空间，其迹为零，故 $\operatorname{Tr}(dP\wedge dP\wedge dP)=0$。又 $dP\wedge dP\wedge dP$ 是 $d(P\,dP\wedge dP)=dP\wedge dP\wedge dP$，所以 $dF_P=0$。十五个分量相加，$dF=0$。
- **曲率与径向条件.** $K_a$ 凸，关于 $0$ 星形。对闭二形式 $F$，Poincaré 引理的同伦公式给出一形式 $\mathcal L(\eta)=\int_0^1t\,F_{tA}(A,\eta)\,dt$，满足 $d\mathcal L=F$。由 $F(A,A)=0$，$\mathcal L_A(A)=0$。这与射线定相一文定理 2 第 3 条的公式相同。
- **适用于命题 5.** 命题 5 第 1 条只用到 $F_P$ 与投影的规范协变性 $P(gU)=gP(U)g^{-1}$；第 2 条只用到上式；第 3 条用到 $d\mathcal L^{\rm rad}=F$ 与常向量场可交换；判据的两个方向只用到第 1–3 条与 $K_a$ 内部连通。投影在 $K_a$ 上的光滑性由引理 1 的隙与 Riesz 积分给出。

**引理 1.**

- **自由核.** 在动量空间（$p_\mu\in2\pi\mathbb Z/L_\mu$；长为一的方向只有 $p_\mu=0$），$D_W(p)=3-\sum_\mu\cos p_\mu\pm i\sum_\mu\gamma_\mu\sin p_\mu$。由 $\gamma$ 互相反交换，

$$
D_W^\dagger D_W=\Bigl(3-\sum_\mu\cos p_\mu\Bigr)^2+\sum_\mu\sin^2p_\mu .
$$

  记 $d_\mu=1-\cos p_\mu\ge0$，则 $3-\sum\cos p_\mu=\sum d_\mu-1$，$\sin^2p_\mu=d_\mu(2-d_\mu)$，

$$
D_W^\dagger D_W-1=\Bigl(\sum_\mu d_\mu\Bigr)^2-\sum_\mu d_\mu^2=2\sum_{\mu<\nu}d_\mu d_\nu\ge0 .
$$

  这对一切取值的 $p_\mu$ 成立，包括 $L_\mu=1,2,3$。又 $\gamma_5D_W\gamma_5=D_W^\dagger$，故 $H(0)^2=D_W^\dagger D_W\ge1$，$|\operatorname{spec}H_q(0)|\ge1$。
- **扰动.** $H_q(A)-H_q(0)=-\gamma_5\sum_\mu B_\mu$，$B_\mu=P^-_\mu\otimes X_\mu+P^+_\mu\otimes X_\mu^\dagger$，$X_\mu=T^U_\mu-T_\mu$。
  - $X_\mu$ 是平移乘以对角阵 $\operatorname{diag}(e^{iqA_{(x,\mu)}}-1)$，故 $\|X_\mu\|\le\max_g|e^{iqA_g}-1|\le|q|\max_g|A_g|\le6a$。
  - $P^\pm_\mu$ 是自旋空间中 $\gamma_\mu$ 的 $\pm1$ 本征子空间上的正交投影，$X_\mu$ 只作用于格点指标。在 $\gamma_\mu$ 的本征基中 $B_\mu=\operatorname{diag}(X_\mu,X_\mu,X_\mu^\dagger,X_\mu^\dagger)$（按自旋分块），所以 $\|B_\mu\|=\|X_\mu\|\le6a$。
  - $\gamma_5$ 酉，故 $\|H_q(A)-H_q(0)\|\le4\cdot6a=1/200$。
- **Weyl 不等式.** 由 Hermite 矩阵本征值的 Weyl 不等式，$|\operatorname{spec}H_q(A)|\ge1-1/200=199/200$。

**定理 2 的证书方法.** 全部比较在 Arb 球算术中完成，精度 $900$ 位（$3^4$ 为 $700$ 位）；Bauer–Fike 一步用 $128$ 位。

1. **链变量与核.** $A_g$ 为精确有理数，$e^{iqA_g}$ 以 Arb 计算。$P_\pm^\mu$ 的矩阵元属于 $\{0,\pm\tfrac12,\pm\tfrac i2\}$，精确输入。$H_q$ 是精确矩阵的 Arb 包络。
2. **谱区间（Bauer–Fike）.** 取浮点本征分解 $W,\Lambda$，把它们的二进制数值当作精确数。以 Arb 计算 $R=HW-W\Lambda$ 与 $E=\|W^*W-1\|$，并验证 $E<1$；谱范数一律以 $\sqrt{\|M\|_1\|M\|_\infty}$ 上界。
   - 于是 $H=W\Lambda W^{-1}+RW^{-1}$，$\|W^{-1}\|\le(1-E)^{-1/2}$，$\kappa(W)\le\sqrt{(1+E)/(1-E)}$。
   - 由 Bauer–Fike 定理，$H$ 的每个本征值 $\mu$ 满足 $\min_i|\mu-\lambda_i|\le\rho=\kappa(W)\|R\|(1-E)^{-1/2}$。
   - 所以 $\min|\lambda_i|-\rho\le|\mu|\le\max|\lambda_i|+\rho$，$H^2$ 的谱含于 $[\mathrm{lo}^2,\mathrm{hi}^2]$。表中的隙下界即 $\mathrm{lo}$。
3. **符号函数的 Legendre 展开.** 取有理 $a\le\mathrm{lo}^2$ 与有理 $s\in(0,1)$，使 $b=a\bigl(\tfrac{1+s}{1-s}\bigr)^2\ge\mathrm{hi}^2$。令 $c=\tfrac{a+b}2$，$d=\tfrac{b-a}2$，$\kappa=c/(1+s^2)$，$u=(t-c)/d$。
   - 由 $b/a=\bigl(\tfrac{1+s}{1-s}\bigr)^2$ 得 $d/c=2s/(1+s^2)$，于是 $t=\kappa(1+s^2+2su)$。
   - 由 Legendre 生成函数 $(1-2xs+s^2)^{-1/2}=\sum_nP_n(x)s^n$（$|x|\le1$）得

$$
t^{-1/2}=\kappa^{-1/2}\sum_{n\ge0}(-s)^nP_n(u),\qquad t\in[a,b].
$$

   - 由 $|P_n|\le1$，截断到 $n\le N$ 的误差一致不超过 $\kappa^{-1/2}s^{N+1}/(1-s)$。
   - $\operatorname{sign}H=H(H^2)^{-1/2}$，$\|H\|\le\sqrt b$，所以 $\|\operatorname{sign}H-Hp_N(H^2)\|\le\sqrt b\,\kappa^{-1/2}s^{N+1}/(1-s)$，$p_N$ 为截断级数。
4. **对角元.** 对标准基向量 $e$，由三项递推 $(n+1)P_{n+1}(X)=(2n+1)XP_n(X)-nP_{n-1}(X)$，在 Arb 中计算 $e^*Hp_N(H^2)e$，$X=(H^2-c)/d$，再把上述截断误差作为球半径加上。$N$ 取使截断误差不超过 $10^{-30}$（$3^4$ 为 $10^{-25}$）的最小者，实际 $N\le247$（$3^4$ 为 $N\le155$）。
   - 由 $\operatorname{tr}_{\rm spin}P_+(x,x)=2+\tfrac12\operatorname{tr}_{\rm spin}\operatorname{sign}H(x,x)$ 得到所需的迹。
   - 零场值同法计算；零场核与 $q$ 无关。
5. **比较.** 由 $\mathcal A$ 的球的上下端点与表中的十进制数比较；“不含零”即球的下端点大于零或上端点小于零。

**推论.** 由定理 2，在三个格上，$\mathcal A(x)(A_*)\ne0$ 对某个 $x$ 成立，且 $A_*\in K_a$。由引理 0，射线定相一文命题 5 对 $\mathcal L^{\rm rad}$ 成立：在 $K_a$ 的内部，$\delta_{\nabla\omega}\mathcal L^{\rm rad}_\eta=0$ 对一切 $\omega,\eta$ 成立，当且仅当 $\mathcal A\equiv0$ 于 $K_a$。命题 5 所需的处处有隙由引理 1 给出。所以前者不成立。

**命题 2.** 由引理 1，在 $K_a$ 上每个 $H_q$ 有隙，$t\mapsto P_+^{(q)}(tA)$ 连续，故秩 $\operatorname{Tr}P_+^{(q)}(A)$ 为常数。所以

$$
\sum_x\mathcal A(x)(A)=3\sum_qw_qq\bigl[\operatorname{Tr}P_+^{(q)}(A)-\operatorname{Tr}P_+^{(q)}(0)\bigr]=0 .
$$

在 $1^4$ 上只有一个格点，故 $\mathcal A\equiv0$。

**定理 3.**

- **第 1 条.** 超荷 $q$ 表示的链变量是 $e^{iqA_g}=e^{i(qA)_g}$，所以 $H_q(A)=H_1(qA)$。由引理 1，$H_1$ 在 $|A_g|\le6a$ 上有隙，Riesz 投影对链变量解析，故 $f$ 在 $0$ 的邻域内是收敛幂级数。对 $|q|\le6$，$qA$ 仍在该邻域内，比较齐次分量得 $f(qA)=\sum_kq^kf_k(A)$。代入 $\mathcal A$ 的定义，零场项消去 $k=0$，即得所写展开。$M_m$ 由整数直接计算。
- **第 2 条.** 取 $C=\gamma_2\gamma_4$。在本文的手征基中，$C$ 是元素为 $0,\pm1$ 的实正交矩阵，满足 $C\gamma_\mu C^{-1}=-\gamma_\mu^T$（$\mu=1,\dots,4$），且与 $\gamma_5$ 可交换；脚本以精确整数矩阵核对这三条。
  - 于是 $C(P^\mp_\mu)^TC^{-1}=P^\pm_\mu$。
  - 又 $(T^U_\mu)^T=(T^{U^*}_\mu)^\dagger$，$\bigl((T^U_\mu)^\dagger\bigr)^T=T^{U^*}_\mu$，转置同时作用于自旋与格点指标。
  - 所以 $C\,D_W(U)^T\,C^{-1}=D_W(U^*)$。
  - 由 $\gamma_5^T=\gamma_5$，$D_W(U)^T=H(U)^T\gamma_5$，故 $H(U^*)=V\,H(U)^T\,V^{-1}$，$V=C\gamma_5$ 是只作用于自旋的酉矩阵。
  - 于是 $\operatorname{sign}H(U^*)=V(\operatorname{sign}H(U))^TV^{-1}$，其 $(x,x)$ 自旋块为 $V_s\bigl(\operatorname{sign}H(U)_{xx}\bigr)^TV_s^{-1}$，自旋迹不变。
  - $U^*$ 对应 $-A$，所以 $f(-A)=f(A)$，奇数阶 $f_k$ 为零。
- **第 3 条.** $3M_5=12960$。

**命题 4.** $L=1$ 时 $(\nabla\omega)_{(x,\mu)}=\omega(x+\hat\mu)-\omega(x)=0$。坐标置换的等价性与射线定相一文定理 2 协变性段的论证相同：超立方转动把 $L_\mu\times L_\nu$ 格映为 $L_\nu\times L_\mu$ 格，自旋转动与 $\gamma_5$ 可交换，所以 $\operatorname{tr}P_+(x,x)$ 随之协变。其余由推论、命题 2 与浮点读数给出。

## 关键数字

- **场.** $A_*=T/5000$，$1/5000<1/4800$。方格角上界为 $1/2500$（$2^4$、$2\times2\times1\times1$）与 $1/1250$（$3^4$）。$6\cdot\tfrac1{1250}\le\tfrac1{200}<\tfrac1{30}$。
- **隙.** 引理 1 给出 $199/200$，即 $1-4\cdot6/4800$。五种超荷中最小的认证 $\mathrm{lo}$：$2^4$ 为 $0.999100$，$2\times2\times1\times1$ 为 $0.998330$，$3^4$ 为 $0.999086$。
- **手征密度.**
  - $2^4$ 格点 $0$：$(2.0362339,2.0362340)\times10^{-11}$。
  - $2^4$ 十六格点：$\min_x|\mathcal A(x)|>2.44\times10^{-12}$。
  - $2\times2\times1\times1$ 格点 $0$：$(-2.6783766,-2.6783765)\times10^{-11}$。
  - $3^4$ 格点 $30$：$(-9.0757235,-9.0757234)\times10^{-13}$。
- **半场比值（$2^4$）.** $\mathcal A(0)(A_*/2)\in(1.2726669,1.2726670)\times10^{-12}$，$\mathcal A(0)(A_*)/\mathcal A(0)(A_*/2)\in(15.9997,15.9998)$。这与定理 3 的四阶首项相符。
- **超荷矩.** $M_1,\dots,M_5=0,-72,0,-1944,4320$；首项系数 $3M_5=12960$。
- **一致性核对（$2^4$）.**
  - $\sum_x\mathcal A(x)(A_*)$ 的球含零。
  - 对 $\omega(i)=((i\bmod3)-1)/5000$ 作规范变换后，$\mathcal A(0)$ 的球与原值重叠。变换后的场不在 $K_a$ 内，但方格角不变。
  - $f(A_*)(0)$ 与 $f(-A_*)(0)$ 的球重叠。

## 与已有结果的关系

- **射线定相一文的命题 6 之后的数值条目.** 原句为：“由于本征分解是浮点计算，结论只是数值：$\mathcal A$ 在四阶非零，指向射线定相不在容许参考类内，未证明。”本文推进了三点：
  - 在 $2^4$、$2\times2\times1\times1$、$3^4$ 上，把“$\mathcal A$ 非零”升级为 Arb 区间证书。所用的场都在 $K_a$ 内部，因而不需要命题 6 与射线上的隙。
  - 把“四阶”升级为定理 3：奇数阶与二阶项恒为零，首项为 $12960\,f_4$。
  - 由引理 0，在这三个格上把“射线定相”定义为不依赖参考测度的连接。
- **该条目的数值.** 其中的三个区间是每个标准块的值，即 $\mathcal A/3$。按该文定义，$\max_x|\mathcal A|$ 约落在 $(6.30,6.33)\times10^{-6}$、$(3.96,3.99)\times10^{-7}$、$(2.478,2.484)\times10^{-8}$ 内；非零性与相邻比值不变。本文的数字按该文定义计入因子 $3$。
- **该文“$3^4$ 是否属于参考测度的体积定义域，也未判定”.** 本文核对了定理 5.4 的陈述：其中只有 $L/\varrho\ge n$，未给出数值，所以该句对 $3^4$、$2^4$ 仍成立。
- **射线定相一文定理 4 的指数局域界.** 它是 $\mathcal L^{\rm rad}$ 这个连接本身的性质，不受影响。本文的结论是：在这三个格上，射线定相不是容许参考测度的候选；在定义域内的体积上它是否是候选，仍未判定。
- **[外尔测度流的曲率部分显式指数局域而对称部分可随参考选择任意放大且衰减任意变慢](外尔测度流的曲率部分显式指数局域而对称部分可随参考选择任意放大且衰减任意变慢.md)的“对所选构造的显式常数仍待给出”.** 不受影响。
- **隙的依据.** 射线定相一文命题 5 原以 Hernández–Jansen–Lüscher 的下界为隙的前提。本文的引理 1 在 $K_a$ 上给出初等的一致隙 $199/200$，推论不依赖该文献结果。
- **复用.** 核 $H_q$、权 $w_q$、$\gamma$ 基与射线定相一文脚本 `ch39-measure-symmetric-part-ray-gauge.py` 中的 `numeric_anomaly_density` 相同；命题 5 的判据直接引用。

## 证据层级

- **书面证明.**
  - 引理 0、引理 1。
  - 推论：由引理 0、引理 1、定理 2 与射线定相一文命题 5 得出。
  - 命题 2。
  - 定理 3 全部三条。
  - 命题 4 的 $L=1$ 部分与坐标置换等价性。
  - 证书方法中的 Bauer–Fike、Legendre 截断误差与 $\operatorname{tr}P_+$ 恒等式。
- **精确有理证书.**
  - 场在 $K_a$ 内部，方格角上界，$1-4\cdot6a=199/200$。
  - 超荷矩 $M_1,\dots,M_5$。
  - $C$ 的三条矩阵性质（精确整数矩阵）。
- **Arb 区间证书.**
  - 定理 2 的全部区间与隙下界，以及隙下界大于 $199/200$。
  - $2^4$ 上十六个格点非零。
  - 半场的值与比值。
- **一致性核对.** $\sum_x\mathcal A$ 的球含零（等式由命题 2 书面证明）；规范变换前后的球重叠（等式由投影的规范协变性书面证明）；$f(\pm A_*)$ 的球重叠（等式由定理 3 第 2 条书面证明）。
- **文献陈述.** Lüscher 定理 5.4 的定义域表述，按 arXiv 版 PDF 的文本提取逐字引用。
- **仅数值.** $2\times1\times1\times1$、$3\times1\times1\times1$、$4\times1\times1\times1$ 上各一个随机场的浮点读数小于 $10^{-13}$。

## 复算

    python3 experiments/ch39-omega-research/replay/ch39-ray-gauge-anomaly-density-certificate.py
    python3 experiments/ch39-omega-research/replay/ch39-ray-gauge-anomaly-density-certificate-three-lattice.py

- 前者给出 $2^4$ 与 $2\times2\times1\times1$ 的证书、一致性核对、半场比值与 $C$ 的精确核对。
- 后者复用前者的函数，给出 $3^4$ 的证书。
- 单核 CPU 时间分别约为一分半与三分钟。

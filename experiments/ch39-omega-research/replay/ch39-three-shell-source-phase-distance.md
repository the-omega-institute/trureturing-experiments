# 截止 52 的同号共门与负号解的源算子距离

固定第五壳源 $S$、四轮 $12,32,52,52$ 的十二条定向弦和六个固定角。令 $U_*$ 取三壳负号解证书的六维活动角盒中的任意一点，最后一轮用 $-G_{52}$。另取任意三个二阶酉门，前两轮与第三轮均可重新选择，第四轮使用第三轮门乘单位相位 $z$，记所得输运为 $V_z$。有
$$
\|U_* S U_*^*-V_z S V_z^*\|_{\rm op}\ \ge\ 6\times10^{-6}|1+z|.
$$
在同号端点 $z=1$ 距离严格大于 $1.2\times10^{-5}$。这个界包括取得三壳六个目标的负号解。其他源算子能否在同号端点给出同样六个读数，仍然开放。[负号角盒](ch39-three-shell-final-sign.md) · [源算子距离的区间复算](ch39-three-shell-source-phase-distance.py)

**从两条重线读出末轮门。** 第五壳源有两条不同的单重重本征线，分别取钉版给出的单位向量 $x_1,x_2$。记 $y_j=U_*x_j$。每条射线的中心、第三轮叶、第四轮叶依次为 $(6,9,13),(7,15,10),(8,16,17)$。在两条本征线上依次列出三条射线的读数，形成六行三列矩阵，三列依次为负号输出的中心 $H$、早叶 $B$、晚叶 $C$。从第三轮前到最终输出，同一个三阶酉矩阵作用于每颗星。故第三轮前后的两套三列矩阵有相同的奇异值。脚本在第三轮前构造矩阵 $M$，取中点矩阵的左逆 $L$，在整个负号角盒上验证
$$
s=\sigma_{\min}[H\ B\ C]\ \ge\ (1-\|I-LM\|_F)/\|L\|_F>0.0061.
$$

在六个叶子输入坐标上倒推候选输运。若候选第三轮门为 $\begin{pmatrix}a&b\\c&d\end{pmatrix}$，从 $H,B,C$ 反推的早叶与晚叶输入分别是
$$
\bar d B+\bar z\bar b(\bar aH+\bar c C),\qquad
\bar z(\bar bH+\bar d C).
$$
前两轮都固定这些输入叶。对负号原门，记其右下元为 $d_*$。减去候选反推时，早叶差在 $B$ 上的系数是 $\bar d_*-\bar d$，晚叶差在 $C$ 上的是 $-\bar d_*-\bar z\bar d$。将其余系数也计入列范数只会增加下界。以 $E^2=\sum_{j=1}^2\|V_z x_j-y_j\|^2$ 表示两向量在原相位下的总误差平方，由酉不变性及限制到各叶的估计，
$$
E^2\ge s^2\left(|\bar d_*-\bar d|^2+|-\bar d_*-\bar z\bar d|^2\right)
\ge \tfrac12 s^2 |d_*|^2 |1+z|^2 .
$$
最后一步是对一个任意复数 $\bar d$ 取平方和的最小值。

**从本征线回到源算子。** 两个输运都固定站点 $14$。令 $p_j=|x_j(14)|>0$，$\Delta$ 为两条重线各自到其余谱的间距的共同下界，$\epsilon$ 为两共轭算子的范数距离。对 $y_j$ 应用候选算子减去第 $j$ 个重本征值，残差范数至多 $\epsilon$。谱分解给出该向量到候选本征线的投影余量至多 $\epsilon/\Delta$，故可选单位相位 $\eta_j$ 使 $\|V_zx_j-\eta_j y_j\|\le\sqrt2\epsilon/\Delta$。在站点 $14$ 取分量得到 $|1-\eta_j|p_j\le\sqrt2\epsilon/\Delta$。于是
$$
E\le\sqrt2\,\epsilon\,C_0,\qquad
C_0=\left(\sum_{j=1}^2[(1+1/p_j)/\Delta]^2\right)^{1/2}.
$$
钉版输入给出 $\Delta>0.3259$、$C_0<272.20$。将同一角盒上的 $s|d_*|/C_0$ 用区间运算直接包围，所得下端大于 $0.000012$。与叶坐标的下界合并即得结论。当 $z=-1$ 且取原门时，可以完全重合。

脚本使用 Python 3.12、python-flint 0.9.0、NumPy、SciPy 和 SymPy。唯一参数为输入提交 9e2b8abfe862acb15667d5a0154731a67e491a95 的 experiments/ch39-omega-research/replay 目录。[第五壳源](https://github.com/the-omega-institute/trureturing-experiments/blob/9e2b8abfe862acb15667d5a0154731a67e491a95/experiments/ch39-omega-research/replay/ch39-fifth-full-positive-joint.py) · [定向弦序与门](https://github.com/the-omega-institute/trureturing-experiments/blob/9e2b8abfe862acb15667d5a0154731a67e491a95/experiments/ch39-omega-research/replay/ch39-two-shell-ordered-bond-gates.py)

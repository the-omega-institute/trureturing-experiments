# 四壳共门的第十一壳留出检验

在平方弦长截止 52、末轮负号的四轮电路中，用同一个第五壳源精确校准第七至第十壳。保持原来的四个固定角，八个活动角在原认证盒中任意变化，在第十一壳的原生第三枢轴得到
$$
0.37799<K_{11}<0.37885,\qquad 0.06957<W_{11}<0.06970 .
$$
两项读数都与事前登记的方向相反。牛顿读数低于校准值 1，真空读数高于校准值 $10^{-122}$。这也覆盖前三壳配平、第十壳目标在原认证矩形中变动的每个盒内解。[四壳校准与目标矩形](ch39-four-shell-final-sign.md) · [留出检验的预登记](https://github.com/the-omega-institute/trureturing-experiments/issues/296) · [区间复算](ch39-four-shell-eleventh-holdout.py)

**留出读数。** 第七壳对应正零点序号 $2\,306\,888\,100\,092$，第十一壳对应 $2\,306\,888\,100\,096$。在后者重新生成荷六壳算子的谱基、第三枢轴的带号权和场表。校准电路的最后一轮在六个出射站点乘以负号，将这些符号乘入谱基伴随矩阵的相应列，再投影输运后的两条重本征线及其导数。若角盒中心为 $c$，半宽为 $r_j$，在中心按热积分计算残差 $F(c)$，在全盒 $B$ 包围雅可比。对每项，在 $F_i(c)$ 两边加减 $\sum_j r_j\sup_B|\partial_j F_i|$，再分别加回目标 $1$ 和 $10^{-122}$。

**真空读数还排除了算子邻域。** 若 $A$ 为角盒中的任一输运后源算子，$A'$ 为同一源的任意酉共轭，且 $\|A'-A\|_{\rm op}\le 0.009$，则 $A'$ 的第十一壳真空读数仍大于 $10^{-122}$。因此即使允许改变输运拓扑，任何达到五壳十式的同谱算子也必须离开这个邻域。

具体地，第十一壳谱基第 $i$ 个单位向量上的源读数为 Rayleigh 商 $x_i$。算子变化至多 $\rho=0.009$ 时，该读数变化至多 $\rho$。把角盒中 $x_i$ 的区间向两端各扩张 $\rho$，再与源算子的谱区间相交。在所得正区间上按端点包围每个因子 $\exp(-m_f/x_i^2)$。若 $w_i$ 是该枢轴的权，$g_f,m_f$ 是钉版场表的带号系数与非负质量参数，则真空项对 $x_i$ 的导数为
$$
-\frac{w_i x_i^3}{2\pi}\sum_f g_f\exp(-m_f/x_i^2).
$$
逐场加总后才取该导数绝对值的区间上界。将各站上界求和并乘 $\rho$，所得变化小于 $0.0631$，小于原角盒与真空目标的间距。角盒外能否在同一四轮电路中另选五壳十式的共门仍开放。

脚本使用 Python 3.12、python-flint 0.9.0、NumPy、SciPy 和 SymPy。唯一参数是输入提交 9e2b8abfe862acb15667d5a0154731a67e491a95 的 experiments/ch39-omega-research/replay 目录。[原生壳谱与第三枢轴](https://github.com/the-omega-institute/trureturing-experiments/blob/9e2b8abfe862acb15667d5a0154731a67e491a95/experiments/ch39-omega-research/replay/ch39-fifth-seventh-shortest-site-gates.py) · [热矩及导数](https://github.com/the-omega-institute/trureturing-experiments/blob/9e2b8abfe862acb15667d5a0154731a67e491a95/experiments/ch39-omega-research/replay/ch39-shortest-common-complex-gate.py)

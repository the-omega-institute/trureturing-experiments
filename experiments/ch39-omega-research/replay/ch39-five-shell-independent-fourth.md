# 四轮独立短弦门上的五壳真空校准

取定第五壳源，在第七至十一壳各自的原生第三枢轴上测量牛顿读数 $K_j$ 和真空读数 $W_j$。仍用截止 $52$ 的十二条定向弦，四轮依次作用在平方弦长 $12,32,52,52$ 上。每轮的三条弦共用一枚 $U(2)$ 门，第四轮可独立选门，整套四枚门在五个壳上共用。对 $t=0.70$ 或 $0.72$，各有一套门同时给出
$$
K_7=K_8=K_9=K_{10}=1,\qquad K_{11}=t,\qquad
W_7=W_8=W_9=W_{10}=W_{11}=10^{-122}.
$$
两个 $t$ 对应不同的门。五个真空目标和四个牛顿目标均精确达到。上述构造给出第十一壳牛顿读数的两个可达到值，原来的 $K_{11}=1$ 仍开放。[角盒数据](ch39-five-shell-independent-fourth.json) · [独立复算](ch39-five-shell-independent-fourth.py) · [三枚门与末轮固定负号的条件结果](ch39-five-shell-conditional-targets.md)

## 末轮多出的可观测自由度

记 $G_1,\ldots,G_4$ 为按实际轮序使用的门。每颗星上，三条共中心的弦依次接入三个不同叶站。第二轮在与星不交的三条弦上作用。若 $D(u)=\mathrm{diag}(u,1)$，在第一轮与第三轮之间、第三轮与第四轮之间给中心引入单位相位 $u,v$，同一个站点输运有以下表示
$$
\bigl(D(u)G_1,\ G_2,\ D(v)G_3D(u)^{-1},\ G_4D(v)^{-1}\bigr).
$$
相邻两轮的中心相位在矩阵乘积中消去。反过来，第一轮叶站的输出行确定 $G_1$ 的第二行，酉性将第一行确定到一个单位相位。对齐并消去第一轮后，第三轮的叶站同样确定 $G_3$ 到第二个中心相位。余下的 $G_4$ 和不交的 $G_2$ 可直接读出。这样也覆盖对角门与交换门。

因此，一套输运能写成末轮等于 $\lambda G_3$（$|\lambda|=1$）的充要条件，是两枚门的右下角模相同。必要性还可直接从叶站 $9,13$ 读出，它们在完整输运矩阵中的对角元分别为 $(G_3)_{22},(G_4)_{22}$。充分性来自两个允许的中心相位。对固定 $\lambda$，若两门右下角满足 $(G_4)_{22}=\lambda (G_3)_{22}$，就可将 $\lambda G_3$ 第一行与第一列分别乘单位相位，使两个非对角元匹配 $G_4$；当非对角元非零时，酉性确定剩下的左上元（右下元为零时左上元也为零）。对角门用一个相位即可。右下元模相同即可选择满足条件的 $\lambda$。

在两张解盒上，右下元模差的一半均严格大于 $0.01$。对任何第三、四轮仅差相位的完整输运，这两个叶站的对角元模相同。由反三角不等式和算子范数控制单个矩阵元，该输运与解盒内任一输运的算子距离均大于 $0.01$。这定量区分了本构造使用的电路；其他输运能否在仅差相位的子族中取得同样的读数仍开放。

## 十项读数的同一组角

五个壳各用本壳的零点、谱基、场表和枢轴权。零点编号依次为 $2\,306\,888\,100\,092$ 至 $2\,306\,888\,100\,096$。同一电路作用在第五壳的两个重列上，投影到各自谱基后用同一组截止谱级计算两项热矩。对每套门固定十六角中的六角，解余下十式。角盒内唯一性只对这些固定角和这个角盒成立。[钉版有序门](https://github.com/the-omega-institute/trureturing-experiments/blob/9e2b8abfe862acb15667d5a0154731a67e491a95/experiments/ch39-omega-research/replay/ch39-two-shell-ordered-bond-gates.py) · [原生谱与第三枢轴权](https://github.com/the-omega-institute/trureturing-experiments/blob/9e2b8abfe862acb15667d5a0154731a67e491a95/experiments/ch39-omega-research/replay/ch39-fifth-seventh-shortest-site-gates.py)

以 $F$ 表示这十项残差，$c$ 为中心，$r_i$ 为半宽，$P$ 为经区间算术核验可逆的常矩阵。以下各量逐行用球算术向外包围
$$
b_i=|(PF(c))_i|,\quad d_i=\sum_j |(I-PDF(c))_{ij}|r_j,\quad
h_i=\sum_{a,b}\sup_{\rm box}|\partial_a\partial_b(PF)_i|\,r_ar_b .
$$
在全十维角盒上界定 $55$ 个无序偏导对，每个混合项按两个有序对计。Taylor 公式给出自映条件 $b_i+d_i+h_i/2<r_i$，导数变化给出压缩条件 $d_i+h_i<r_i$。各行同时满足两条件，故存在一个精确零点。小目标 $10^{-122}$ 用原生球算术计入残差，两项的证据均为角盒内的存在证明。[钉版热核与二阶导数](https://github.com/the-omega-institute/trureturing-experiments/blob/9e2b8abfe862acb15667d5a0154731a67e491a95/experiments/ch39-omega-research/replay/ch39-three-shell-gravity-radius56.py)

复算使用 Python 3.12、python-flint 0.9.0 和输入提交 `9e2b8abfe862acb15667d5a0154731a67e491a95` 的 replay 目录。

```sh
python3 ch39-five-shell-independent-fourth.py /path/to/pinned/experiments/ch39-omega-research/replay
```

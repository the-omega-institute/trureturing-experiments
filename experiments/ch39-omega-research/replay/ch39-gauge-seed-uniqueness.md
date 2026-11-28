# Uniqueness of the sharp two-loop crossings used in the portal fit

The cut-9 clock model uses the sixteen-scalar spectrum, denominator $724$, the fixed center-even mass multiplier $32$, and sharp activation of each center-odd mass. The three inverse gauge couplings and the three logarithmic Yukawa couplings evolve by the [published two-loop gauge and one-loop Yukawa equations](ch39-two-tier-dirac-running.py). For each of the [four fit and ten holdout candidates](ch39-cubic-portal-holdout.py), there is exactly one crossing in
$$
 28.940\le L=\log(\mu/M_Z)\le29.140,\qquad
 0.117\le\alpha_s(M_Z)\le0.119 .
$$
The published [seed enclosure](ch39-chiral-clock-unit-readout.py) proves existence in the smaller rectangle $L=29.040\pm0.014$, $v=1/\alpha_s(M_Z)=1/0.1184\pm0.04$. In particular the enclosed solution also lies in the larger rectangle. The uniqueness argument treats the changes of the beta functions at mass thresholds. Both rectangles stay below every lifted center-even threshold.

## Crossing equations with a movable endpoint

Write $X_i=1/\alpha_i$, $h_j=\log y_j$, $g_i=1/X_i$ and $Y_j=e^{2h_j}$. Between two sharp mass thresholds the equations have the form
$$
 X_i'=-b_i/(2\pi)+f_i,\qquad
 f_i=-\sum_j B_{ij}g_j/(8\pi^2)+\sum_j D_{ij}Y_j/(32\pi^3),\qquad
 h_i'=\frac{\sum_j A_{ij}Y_j-4\pi\sum_j C_{ij}g_j}{16\pi^2}.
$$
Here $2A$ has rows $(9,3,2),(3,9,2),(6,6,5)$. Denote $c_i(L,v)=\int_0^L f_i(t,v)\,dt$. In either rectangle let $b_i^0$ include the thresholds strictly below its lower face. The seed map solves a fixed linear system whose rows are $(b_i^0,2\pi,-2\pi\mathbf1_{i=3})$ for $(L_{\rm out},X_{\rm common},v_{\rm out})$. Denote the matrix by $K$. Thresholds possibly opening inside the rectangle contribute the terms $-\Delta b_i(L-t_e)_+$ on the right. The remaining terms there are $2\pi c_i$ and constants. At a smooth endpoint, the two columns of the derivative are therefore obtained by applying $K^{-1}$ to
$$
 2\pi f_i(L,v)-(b_i(L)-b_i^0),\qquad
 2\pi\left[\frac{\partial X_i(L,v)}{\partial v}-\mathbf1_{i=3}\right].
$$
Discard the middle output row to get the derivative of the map on $(L,v)$. At each threshold the solution is continuous and the Jacobian bounds from the two adjacent open intervals give a Lipschitz bound across the join.

For the $v$ derivative integrate $s'=J s$, $s(0)=(0,0,1,0,0,0)$, alongside the published flow. Here $J$ is the Jacobian of its six velocities in $(X,h)$. Its four blocks have entries $B_{ij}g_j^2/(8\pi^2)$, $D_{ij}Y_j/(16\pi^3)$, $C_{ij}g_j^2/(4\pi)$ and $A_{ij}Y_j/(8\pi^2)$ respectively. The mass thresholds are fixed when $v$ varies. On each step the existing interval integrator supplies enclosures for $g$ and $Y$ on the entire step. Widen every sensitivity coordinate by one quarter of the upper step length and verify that the velocity interval times the step fits strictly inside this padding. The integral-equation enclosure then contains all intermediate sensitivities and its forward endpoint update encloses the true endpoint. Endpoints in the rectangle are tested throughout each overlapping step. Integrating to the upper terminal time covers every possible value of $L$.

In the norm $\max(|\delta L|/0.014,|\delta v|/0.04)$ the largest row sum of the rescaled derivative is below $0.152$ on the original seed and below $0.42$ on the larger rectangle for all fourteen candidates. Thus two distinct fixed points in either rectangle are impossible. This last inference only uses the Lipschitz bound; existence comes from the smaller seed's image enclosure. The cutoff scheme, spectrum and candidate index are fixed separately in each comparison. Scale ranges outside the larger rectangle require their own exclusion test. [Fourteen interval derivative checks](ch39-gauge-seed-uniqueness.py)

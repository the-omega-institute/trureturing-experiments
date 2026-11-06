# Normalizer and centre-index check

`check.py` implements the groups
\(S=\langle a,b\mid a^8=b^2=1,\ b^{-1}ab=a^3\rangle\),
\(P_p=\langle x,y\mid x^{p^2}=y^p=1,\ yxy^{-1}=x^{1+p}\rangle\)
with integer pair multiplication. It enumerates all subgroups of `S`, `P_3`, and `P_5`; for the coprime product `S × P_3` it enumerates all Cartesian products of the two coordinate subgroup lists, as supplied by the coprime direct-product splitting theorem, and computes normalizers and their conjugacy orbits exactly.

Run from this directory:

```text
python3 check.py
```

The stored output is `check.out`; the final line is `ALL_CHECKS_OK`. The checks reproduce `s(S)=5`, `ω(S)=4`, `s(P_p)=ω(P_p)=2`, and for `S × P_3`, `s=10`, `ω=8`, and `|G:Z(G)|=72=8·3²`.

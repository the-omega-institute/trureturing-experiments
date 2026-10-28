# Rehman et al. Turán–Seidel insertion check

This exact SymPy check accompanies the preregistration for Question 1 of Rehman, Yang, Ajmal, Rodrigues, Maden and Cai, *Seidel polynomials of an edge embedded t-partite Turán graphs and Seidel energy analysis for the 4-partite case*, AIMS Mathematics 11(10) (2026), 32514–32537, DOI 10.3934/math.20261279.

`check.py` uses integer matrices and exact polynomial/radical arithmetic only. It verifies the source's balanced `p=2` family (`T(2t,t)`), including the factorisations and the strict negative energy increment for `t=3,…,8`; it checks the independent-twin threshold algebra `2−8/m`; and it checks the proposed sharpness polynomial with `det(4I−M)=3(t−7)`, giving equality at `t=7` and a decrease for every `t≥8`. The final line of the stored output is `ALL_CHECKS_OK`.

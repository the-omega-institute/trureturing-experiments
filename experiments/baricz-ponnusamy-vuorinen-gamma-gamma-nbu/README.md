# Baricz–Ponnusamy–Vuorinen Question 8: is the gamma-gamma distribution NBU?

Companion programs for the negative answer to Question 8 of Á. Baricz, S. Ponnusamy and M. Vuorinen, *Functional inequalities for modified Bessel functions*, Expo. Math. 29(4) (2011) 399–414, arXiv:1009.4814.

    python3 gg-cert.py    # rigorous outward-rounded rational quadrature: S(1) > S(1/2)^2 at a = b = alpha = 1 (standard library)
    python3 gg-scan.py    # mpmath check of the density, the survival closed form and S(2u) - S(u)^2 on a grid

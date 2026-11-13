# Asymptotic spanning tree candidate check

`check.py` uses integer-only arithmetic. It enumerates reduced words in a radius-five ball of the Cayley graph of `Z^2 * Z`, checking inverse-closed generators, preservation of the first free-product syllable away from the identity, and the isometric `Z^2` factor distances. It then exhausts all spanning trees of the 2x2 and 3x3 square grids, measuring the maximum ambient distance of the tree path over every grid edge. The 3x3 minimum edge stretch is 3.

The finite grid enumeration is a bounded check supporting the one-ended-factor obstruction; the proposed infinite refutation is the mathematical argument in the preregistration.

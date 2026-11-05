# Neumann partition separation

This report checks the positive-weight complete-graph family proposed for the open comparison question in Section 5.1 of Hofmann, Kennedy, Mugnolo and Plümer, arXiv:2608.19962v1. `check.py` uses SymPy exact rational algebra only. It expands the source Neumann form on every two-vertex mean-zero cluster, enumerates all seven unordered exhaustive bipartitions of `K4`, verifies the strict Neumann ordering for every real parameter `t >= 6`, and checks the boundaryless and Dirichlet certificates. The final output line is `ALL_CHECKS_OK`.

# Logarithmic norm counterexamples

Source: M. Hayajneh, S. Hayajneh and F. Kittaneh, *On some classical trace inequalities and a new Hilbert–Schmidt norm inequality*, MIA 21(4) (2018), 1175–1183, DOI [10.7153/mia-2018-21-81](https://doi.org/10.7153/mia-2018-21-81), Conjecture 2 and Theorem 3.

Run `python3 check.py` with Python 3 and SymPy. The expected output is in `output.txt`; its final line is `ALL_CHECKS_OK`.

The program uses symbolic polynomial identities and rational interval arithmetic, without floating-point evidence. It reproduces Theorem 3 on the family A=diag(exp(-1),1), B=[[1,b],[b,1]], including equality at t=1/2. Put u=exp(-t), v=exp(t-1), c=v-u. For 1/2<t<1 and c/[v(1-u)]<b²<1, the trace-norm comparison is reversed, with squared gap 4c(1-b²)>0, while the Hilbert–Schmidt squared gap in the source's direction is 2cb²>0.

Exact Taylor bounds certify t=3/4, b=7/8 uniformly for A_eps=diag(exp(-1+eps),exp(-eps)), 0<=eps<=1/10000: det(L)<-3/8, det(R)<-1/200, and the trace-norm squared gap exceeds 1/4. For positive eps in that interval, both eigenvalues lie strictly inside (exp(-1),1), and B is positive definite. The endpoint eps=0 supplies a concrete counterexample to Conjecture 2.

This is a mathematical check program, not a Lean proof. It does not certify the separate proposed Schatten-p asymptotic extension or a literature-priority claim. No source PDF or copyrighted source text is distributed here.

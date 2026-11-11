# Finite metric doubles for the Met terminal-chain question

`check.py` uses integers and an explicit positive-infinity sentinel. It enumerates every separated symmetric extended metric on up to four labelled points whose positive distances belong to {1, 2, 3, +infinity}. For every nonempty subspace it checks the proposed double, its isometric copies, nonexpansive fold, intersection square, and the resulting pullback distance identity for categorical powers of arity zero, one and two. A separate shortest-path calculation checks that the proposed distances admit no shorter paths.

The fidelity check reproduces the common-stage distance formula for directed unions of subspace embeddings stated in Adámek–Milius–Moss, *Terminal Coalgebras in Countably Many Steps*, arXiv:2303.11071v4, proof of Lemma `L:emb`, item 1 (the corresponding finite-subspace factorization lemma in CALCO 2025, §6).

Run with Python 3.12 or newer:

    python3 docs/reports/adamek_milius_moss_met_omega_plus_omega/check.py

`output.txt` contains the deterministic output. These are finite sanity checks of the proposed reduction. They do not prove preservation for arbitrary endofunctors, preservation of directed colimits, or convergence at omega + omega. The proposed general proof must establish those steps separately.

Source: Jiří Adámek, Stefan Milius and Lawrence S. Moss, *Terminal Coalgebras for Finitary Functors*, CALCO 2025, LIPIcs 342, article 3, DOI [10.4230/LIPIcs.CALCO.2025.3](https://doi.org/10.4230/LIPIcs.CALCO.2025.3), §7, second open question. The combined arXiv v4 retains the theorem with the extra isometric-embedding hypothesis but does not repeat this open question.

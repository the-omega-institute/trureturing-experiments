# Exact HS4 exchange readout and observable closure

The program `verify_exchange_readout.py` supports the ordinary mathematical [HS4-AURIC network companion](https://github.com/the-omega-institute/trureturing/blob/4e51f7b21bedd7586628e541e24e7d6caa13c74d/docs/develop/theory/HIGUCHI_SUDBERY_AURIC_NETWORK_CLOSURE.md). It requires Python and SymPy. The reference environment is Python 3.13.5 and SymPy 1.14.0. It makes no network calls and writes JSON to standard output.

```sh
python verify_exchange_readout.py > actual.json
python -O verify_exchange_readout.py > optimized.json
cmp actual.json optimized.json
python -m py_compile verify_exchange_readout.py
```

The reference output is `exchange_readout_result.json`: 190 exact finite/symbolic checks and three negative controls. The two executions above produced identical JSON. The source uses explicit exceptions rather than assertions, so optimization does not remove the checks.

The exact source SHA-256 is `ab5e4b295dbf955702cceaefb3c7da4f42976170fd90ee585cd39e9fb2459a2c`; its Git blob is `39a9e9f4e4d2a8b9978479b95ef9f5b44d990587`. The uploaded source blob matches the executed local source.

Coverage includes every weight-at-most-two four-qubit Pauli operator and the resulting restricted observation rank; full-space exchange preservation, cubic identity and unitary pulse; all-angle hidden-coordinate response; finite real-projector histories; the FIB symmetric product defect and the logical alternating commutator defect; Fisher and pulse-sign erasure identities; all six HS pair marginals; tetrahedral separable decomposition and balanced-sign moments of the attaining approximation; all 27 triple tensor entries; and all eight probabilities of the independent-source triangle parity construction.

The fully separable distance bound and the arbitrary-length real-instrument no-go are proved by universal arguments in the companion. Finite histories in the checker only test instances. A symbolic or finite check is not a new Lean proof or an independent peer review. This program does not recompile the inherited HS4 proof, certify a source-network model for every measurement setting, optimize arbitrary control menus, or establish a laboratory implementation. Signed simultaneous exchange terms are an explicit control assumption. Earlier experimental programs and their reference outputs are separate preserved evidence.

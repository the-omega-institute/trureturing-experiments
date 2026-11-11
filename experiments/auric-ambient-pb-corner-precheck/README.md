# Independent matrix and symbolic checks

With Python 3.12 or newer, install `requirements.txt` and run:

```sh
python3 experiments/auric-ambient-pb-corner-precheck/independent_precheck.py
```

The program reconstructs Gaussian-integer reflector numerators, adjoints,
full Hermitian input coordinates and every matrix coefficient without importing
or executing the producer. It compares all four arrays in
`../auric-ambient-pb-corner-research/`, checks the oriented determinant separately
by integer modular elimination, and checks full-product ambient and symbolic
strict-center identities. `check-results.json` is its canonical result dataset.
Run the producer first if regenerating its arrays.

Frame and basis definitions come from the owning main theory volume §§48–50;
full-input and ambient applications are in §§53–54. The source pin in the result
identifies the mathematical definitions, not a required local Git object.
All files resolve relative to this checkout, including when run from another
working directory. No second checkout or diagnostic logs are required.

This is executable calculation evidence. It does not certify Lean, physical
source availability, all native contracts, literature priority, or final
independent delivery review. Reuse of the same arithmetic libraries is explicit.

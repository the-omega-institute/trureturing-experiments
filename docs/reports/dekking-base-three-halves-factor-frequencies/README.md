# Dekking base-3/2 factor-frequency checks

`check.py` implements the canonical base-³/² digit extraction using only integers and `fractions.Fraction`. It reproduces Dekking's table (`SQ(7)=211` and the displayed A357448 prefix), checks the defining value identity and the 2-to-3 recurrence on 100,000 exact cases, and checks the two conjectured pair frequencies against a one-million-term prefix with exact rational error bounds. These are finite fidelity and sanity checks; they do not prove a limiting frequency.

Run with:

```sh
python3 check.py
```

The recorded output is in `check.out` and ends with `ALL_CHECKS_OK`.

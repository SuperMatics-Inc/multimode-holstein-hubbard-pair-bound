# Finite replay

This directory contains the deterministic finite checks accompanying the preprint. Install Python 3.12 and `mpmath` 1.4.1, then run from the repository root:

```sh
python -m pip install mpmath==1.4.1
python -B replay/replay.py
python -B replay/replay.py --all-mutants
```

The first command checks the admitted finite witnesses and constants. The second verifies rejection of all 20 guarded mutants. `SOURCE-PINS.json` records SHA-256 digests for the reviewed proof-consumer source files preserved under `ADIABATIC/`.

These finite checks do not formalize the infinite-dimensional oscillator trace identity, the infinite-volume cluster expansion, or the canonical pressure-transfer argument. The mathematical claims and their hypotheses are stated in `paper/main.tex`.

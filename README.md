# A frequency-uniform pair bound for a multi-mode Holstein-Hubbard model

Davis Rens  
SuperMatics Inc.  
[ORCID 0000-0001-5839-5197](https://orcid.org/0000-0001-5839-5197)  
Version 1.0.0

Paper DOI: <https://doi.org/10.5281/zenodo.24384461>  
Software supplement DOI: <https://doi.org/10.5281/zenodo.24384463>  
Source: <https://github.com/SuperMatics-Inc/multimode-holstein-hubbard-pair-bound>

This repository accompanies the preprint *A frequency-uniform pair bound for a multi-mode Holstein-Hubbard model*. The manuscript source and bibliography are in `paper/`. The deterministic finite replay is `replay/replay.py`; its reviewed proof-consumer source files are preserved as byte-identical copies under `replay/ADIABATIC/` and pinned by `replay/SOURCE-PINS.json`. `MANIFEST.sha256` pins every file in this repository except the manifest itself. `ADMITTED.sha256` records the original raw-byte inputs used to prepare the release.

## Mathematical scope

For spinful electrons on square tori of side at least four, with finitely many independent local optical modes of positive frequency, define `u = sum_r F_r^2/K_r` and `U_eff = u - U_c`. When `U_eff > 0` and `beta t^2/U_eff <= 2^-27`, the paper bounds the absolute equal-time on-site pair correlation at distinct sites by `(240/59)(5/64)^d`, uniformly in every positive mode frequency and the number of modes. The associated canonical pair-temperature bound applies to translation-invariant thermodynamic limits at interior densities, not as a uniform finite-volume canonical estimate.

A nonlocal reduction is proved only when the mass-weighted Gram matrix is exactly diagonal for a common-frequency group, with the condition satisfied separately in each equal-frequency group. Static cancellation between unequal-frequency groups is not enough. The counterexamples show failure of one factorization method, not a no-go theorem for generic nonlocal models. A frequency-uniform bound for generic summable nonlocal forces without an exact orthogonal reduction remains unresolved. These are mathematical results for the stated model, not material predictions, synthesis results, or measurements.

## Replay

Use Python 3.12, `mpmath` 1.4.1, and `numpy` 2.5.5. From the repository root:

```sh
python -m pip install mpmath==1.4.1 numpy==2.5.5
python -B replay/replay.py
python -B replay/replay.py --all-mutants
```

The replay uses exact rational arithmetic and 70-digit outward intervals for finite checks. The second command verifies rejection of all 20 guarded mutants. These finite checks do not formalize the infinite-dimensional oscillator trace identity, the infinite-volume cluster expansion, or the canonical pressure-transfer argument. The analytic arguments and their hypotheses are in the paper.

## Disclosure and review

AI agents generated the proofs, code, and drafts under D. Rens' direction. Independent reviews were performed by separate AI agents. This preprint has not been peer reviewed.

## Rights

The manuscript text and bibliography in `paper/` are distributed under CC-BY-4.0 through the paper record. The replay code and this repository's software materials are distributed under Apache-2.0. The paper's license does not change the software license.

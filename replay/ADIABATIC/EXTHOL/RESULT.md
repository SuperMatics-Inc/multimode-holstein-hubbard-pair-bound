# EXTHOL result

Date: 2026-11-19. Assigned files only: `ADIABATIC/EXTHOL/RIGOR.md`, `replay.py`, `RESULT.md`. Workspace is not a Git repository; no branch or commit was created. No frozen lane or paper file was edited. No publication, deployment, repository visibility change, email, cloud job or Git operation was performed.

## 1. Three separate outcomes

### 1.1 On-site multi-mode theorem

**PROVED**, with finite constants **CERTIFIED** in the replay. On a square torus of side at least four, take nearest-neighbor spinful hopping t, finitely many independent on-site optical modes with real F_r, K_r>0 and Omega_r>0, an on-site repulsion U_c>=0, and arbitrary real centered chemical potential. Set

$$u=\sum_rF_r^2/K_r,\quad U_{\rm eff}=u-U_c>0,\quad G=U_{\rm eff}/2.$$

There is no restriction on any ratio of mode frequencies. If beta t^2/U_eff<=2^{-27}, then the distinct-site on-site electronic pair correlation obeys

$$|\langle\Delta_x^\dagger\Delta_y\rangle_{\rm GC}|\le\frac{240}{59}\left(\frac5{64}\right)^{d(x,y)}.$$

The absolute finite-volume pair-matrix row sum is at most 512579/205379. With the reviewed real-source canonical transfer and the same on-site pair-order definition as STAGE7R/HUBHOL, at every density rho in (0,2),

$$k_BT_{\rm pair}\le 2^{27}\frac{t^2}{U_{\rm eff}}=2^{26}\frac{t^2}{G}.$$

Here G is the atomic cost per odd site, attained at zero chemical potential in the local atomic model. It is not mislabeled as the energy cost of two odd sites, which is 2G in that case. The constants are independent of the finite number of branches and of every positive frequency. An infinite branch set, a zero frequency itself, and arbitrary nonlocal pair wavefunctions are not included. See `RIGOR.md` sections 1.3 and 2.

### 1.2 Nonlocal orthogonal-kernel reduction

**PROVED**, with an exact finite orthogonality example **CERTIFIED**. At common frequency Omega>0, suppose on every finite torus

$$\sum_{j,r}\frac{f_r(i-j)f_r(k-j)}{K_r}=u\delta_{ik},\qquad u>U_c.$$

An oscillator-only orthogonal rotation converts the coupled coordinates into independent on-site bright modes; all dark oscillators decouple. The theorem and constants in section 1.1 therefore apply unchanged. The replay checks the four-cycle kernel (1/2,1/2,-1/2,1/2), K=1, whose cyclic Gram matrix is the identity, embedded along rows of the square torus.

For several frequency groups, it is sufficient that each equal-frequency group's Gram matrix separately be u_g delta_ik; then u=sum_g u_g. Cancellation of only the total static off-diagonal Gram entries across unequal frequencies is insufficient for this reduction. The replay certifies a nonzero dynamic cross exponent for two oppositely signed kernels whose total static nearest-neighbor interaction cancels. See `RIGOR.md` section 3.

### 1.3 General nonlocal coupling: a method obstruction, not a model no-go

**PROVED**, with two finite witnesses **CERTIFIED**. The Lang-Firsov transformation and a global nonnegative Gaussian action remain valid. But the STAGE6 factorization into compatible connected physical-support polymer activities fails: disjoint electronic supports need not have independent atomic weights or oscillator traces. Small fixed off-diagonal row sum does not make that equality true.

Electronic Gram witness: U_00=1, nearest horizontal-neighbor U_01=1/16, and Delta_1=1/8<=U_00/4. It is realized exactly by f_1=(1/4,1/4), K_1=1, and f_2=(1,0), K_2=8/7. At beta=1 and t=1/8, two disjoint hopping edges violate Z_02 Z_empty=Z_0 Z_2 in the exact finite electronic reference.

Finite-frequency single-mode witness: f=(1,1/16,0,0), K=1, Omega=6/5, beta=1, with U_00=257/256 and Delta_1=1/8<=U_00/4. Closed hop-and-return trajectories on two disjoint edges at times (1/5,3/5) and (2/5,4/5) have P_02>P_0 P_2.

Both examples embed along rows of a square lattice or four-site-side square torus. They are not no-go theorems for a generic extended Holstein model, another choice of polymer labels, or another expansion. A general Omega-uniform pair-temperature theorem under summability or diagonal dominance alone is **OPEN** in this lane. See `RIGOR.md` section 4.

## 2. Claim ledger

| Claim | Label | RIGOR section | Scope/evidence |
| :-- | :-- | :-- | :-- |
| Summable Gram kernel and correct torus folding | PROVED; CERTIFIED finite witness | 1.1 | l1 convolution estimate and exact rational folding |
| Multi-mode Lang-Firsov atomic energy and centering | PROVED; CERTIFIED finite algebra | 1.1, 1.2 | correct U_eff and chemical-potential shifts |
| Nonlocal atomic parity bound E>=E_even+G N_odd | PROVED; CERTIFIED finite samples | 1.3 | independent random rounding; four complete 81-configuration samples |
| Positive chronological on-site phonon factor, independent supports | PROVED | 2.2 | reviewed infinite-oscillator identity branch by branch; no cutoff |
| Arbitrary-hopping edge SOS and 5/4 cap | PROVED; CERTIFIED finite algebra | 2.3 | exact CAR/SOS and polynomial identity; no global sign gauge |
| Physical cluster constants and pair row sum | PROVED; CERTIFIED | 2.1, 2.4 | q=5/64, A=240/59, physical budget=75/236 |
| Canonical source constants | PROVED; CERTIFIED finite constants | 2.4 | source disk tau/2^14, q_src=45/256, budget=675/844 |
| Canonical on-site pair-temperature conclusion | PROVED | 2.1, 2.4 | reviewed real-source Legendre transfer, every rho in (0,2) |
| Common-frequency orthogonal nonlocal reduction | PROVED; CERTIFIED finite example | 3.1 | oscillator canonical rotation; exact cyclic Gram identity |
| Separate equal-frequency-group reductions | PROVED | 3.2 | rotation within each degenerate group only |
| Static cancellation with unequal frequencies need not cancel dynamics | PROVED; CERTIFIED witness | 3.2 | nonzero outward interval for the cross exponent |
| Failure of STAGE6 atomic polymer factorization under dominance | PROVED; CERTIFIED finite witness | 4.1, 4.2 | actual spinful CAR traces, not just a Gibbs-measure assertion |
| Failure of finite-frequency oscillator support factorization | PROVED; CERTIFIED witness | 4.3 | independently evaluated chronological Gaussian traces |
| Fixed dominance does not ensure beta-independent small Mayer labels | PROVED; CERTIFIED scalar witness | 4.4 | beta=16, U_01=1/16 gives e-1>e^-2/16 |
| Four-site local electronic diagonalization comparison | NUMERICAL | 6 | floating point, not a finite-frequency oscillator simulation |
| Generic nonlocal gap-to-pair theorem | OPEN | 1.3, 4.4 | needs a nonproduct reference or another resummation |

Section 5 screens every transplanted proof ingredient against the STAGE6/STAGE7R reviews. In particular, it uses decreasing-time oscillator order, the SOS constant 2, corrected degree-squared support counts, exact deletion denominators, the genuinely uniform earlier STAGE6 cap narrative, and no false slow-phonon/bare-block identification.

## 3. Observed replay commands and exits

The explicitly assigned command prefix was used on the laptop:

```text
UV_OFFLINE=1 UV_CACHE_DIR=D:/dev/.cache/uv uv run --no-project --with numpy --with scipy --with mpmath --with python-flint python -B
```

Applied to `D:/dev/scratch/ceiling-20261118/ADIABATIC/EXTHOL/replay.py`:

| Final command | Observed exit | Observed result |
| :-- | :-- | :-- |
| `replay.py` | 0 | all exact-rational and outward-interval checks passed |
| `replay.py --all-mutants` | 0 | all 20 named children exited 1 with matching `REJECTED <name>:` prefix |

The final combined invocation exited 0, with `EXTHOL_BASE_EXIT=0` and `EXTHOL_MUTANTS_EXIT=0`; observed wall time was 113.95 seconds. The children were sequential.

All observed mutant exits were 1:

```text
lf-uncentered
ueff-sign
torus-wrap-missing
multimode-sum-wrong
polymer-factor-nonlocal
phonon-factor-nonlocal
mayer-lowtemp
orthogonal-kernel-wrong
charge-gap-double
sos-c1
jw-gauge
cap-54-wrong
scale-too-small
threshold-constant-wrong
physical-degree-wrong
pair-prefactor-wrong
source-disk-wrong
source-degree-wrong
row-shell-wrong
cap-identity-wrong
```

Development history is not hidden: an earlier base run exited 0, but its mutant controller exited 1 because `charge-gap-double` escaped its sampled nonlocal configurations and exited 0. Adding the sharp local zero-chemical-potential witness repaired that guard. The final replay additionally checks the per-odd-site inequality on every configuration, and the repaired mutant now exits 1 with `Odd-sector gap 1/2 < claimed=1`. The final exits above apply to the delivered replay bytes.

### 3.1 Delivered evidence

All interval decisions used outward mpmath.iv enclosures at 70 digits. The final runs reported:

* **CERTIFIED** oscillator difference P_02-P_0 P_2 lies between 0.00308083514745140402766654538067031773665765003191782938119414874823911764593 and 0.0030808351474514040276665453806703177366576500319178293811941487482406027773, strictly positive.
* **CERTIFIED** unequal-frequency static-cancellation cross exponent lies between 0.00208775744739214180885895913819793734548404006343557325488138992760934959069 and 0.00208775744739214180885895913819793734548404006343557325488138992760945606567, strictly positive.
* **CERTIFIED** electronic Z_02 Z_empty-Z_0 Z_2 lies between 0.0126472672543419958292941707749112592315473339737085759951351268759009678931 and 0.0126474669598588226047329425687964559877514294956364779599059718022912130263, strictly positive.
* **NUMERICAL** four-site local electronic Omega=infinity diagonalization reported maximum correlation/bound ratio 7.597e-15 for two on-site Hubbard parameter choices. This is a numerical consistency check, not a certificate of the theorem or a thermodynamic result.

Largest exact interval blocks are 16 by 16; the separate numerical matrix is 256 by 256. No oscillator Hilbert-space truncation or large-memory computation was used. Peak process RSS was not instrumented; no measured RAM figure is claimed.

Other exercised commands: raw-byte `sha256sum` on RIGOR.md and replay.py exited 0; `rmdir D:/t/ExtHolstein` exited 0, with `EXTHOL_CLEANUP_EXIT=0`. The assigned temporary directory was empty before its removal. Builds, repository tests, linters, formatters, TeX/PDF rendering, formal theorem-prover verification and cloud jobs were not run. Only the assigned replay modes were executed.

## 4. File hashes and independent review status

SHA-256 of the delivered source bytes:

* `RIGOR.md`: `c4a38075787a2af8a5aa73ec7aaabf3afc3cc8a72cdae23864ea8e0774965545`
* `replay.py`: `e926e9cc641e6002c0f3acb37cbe1750dbe990575f34c621180357d5402135dd`

The raw-byte SHA-256 of this RESULT.md is supplied in the READY delivery receipt rather than circularly embedded in its own bytes. That receipt supplies all three final file hashes.

**Ready for independent review: yes.** This extension has not yet received its own independent review. Its analytic dependencies have the cited separate reviews; the finite replay modes have the observed passing exits. Independent review should examine the multi-mode chronological trace and source transplant, the exact orthogonal-coordinate reduction, the atomic rounding argument, and the two specific factorization witnesses. The unresolved generic nonlocal theorem is explicitly separated from the completed results and is not a release claim.

## Review and fixes (Main, 20 Nov 2026)
- EXTHOL-REVIEW/REVIEW.md: claim 1 PASS 0.97, claim 2 PASS 0.97, claim 3 FIX 0.98 (one P2 algebra slip, one P3 wording). Replay base exit 0, 20/20 mutants rejected, reviewer mutants cg-225 and g-no-half both rejected.
- Main applied both fixes to RIGOR.md. In §4.3 the bracket is g(1/5)-g(2/5), not 3g(1/5)-g(2/5). In §2.3 the "four CAR sectors" wording was clarified. Neither the conclusion P_02>P_0P_2 nor any replay check changes; replay.py is untouched.
- RIGOR.md sha256 after fixes: b92d7b423e64be51b530c468fc43bc3144652f20674558a1e86cce817be6510d (it supersedes c4a38075... in the hash list above).
- EXTHOL-REVIEW/DELTA.md: PASS 0.99; claim 3 upgraded to PASS. All three claims are now reviewed and pass.
- Main decision: EXTHOL becomes a short standalone preprint (ADIABATIC/EXTHOL-PAPER), released like HUBHOL, so the capstone cites a released DOI and does not carry new proofs.

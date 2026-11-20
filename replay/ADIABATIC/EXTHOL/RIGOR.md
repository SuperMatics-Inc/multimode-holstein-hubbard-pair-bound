# EXTHOL: several optical modes and nonlocal factorization

Davis Rens  
SuperMatics Inc.  
ORCID 0000-0001-5839-5197

This is an internal proof report, not a paper body. Claim labels have the meanings: **PROVED** is an analytic argument below or a reviewed argument whose hypotheses are rechecked below; **CERTIFIED** is a finite exact-rational or outward-interval witness in `replay.py`; **NUMERICAL** is floating-point evidence, not proof; **[INFERENCE]** is a suggested implication without proof; **OPEN** is unresolved here. These labels must not be copied into a preprint body.

AI agents generated the proofs, code and drafts under D. Rens' direction, and separate AI agents did the independent reviews of the prior results used here. This extension is ready to be assigned its own independent review. This preprint has not been peer reviewed.

## 0. Outcome and provenance

**PROVED:** the complete attractive Hubbard-Holstein theorem extends to any finite number of independent on-site optical modes, with unchanged constants and no restriction on their positive frequencies. A second exact reduction covers nonlocal kernels whose mass-weighted coupling rows are orthogonal within each common-frequency group. General nonlocal couplings are not covered by either theorem.

**PROVED:** for genuinely nonlocal Gram interactions, the STAGE6 disjoint-support factorization is false. Sections 4.2 and 4.3 give distinct electronic and finite-frequency oscillator witnesses, both satisfying an explicit diagonal-dominance condition. This is a failure of that factorization, not a no-go theorem for the model or for another expansion.

Sources used, all read-only and relative to `D:/dev/scratch/ceiling-20261118/`:

* `ADIABATIC/LINEAR5/STAGE4/RIGOR.md`, sections 4 and 5: spatial polymer deletion and real-source canonical transfer.
* `ADIABATIC/LINEAR5/STAGE6/RIGOR.md`, sections 2, 3, 6 and 7: chronological oscillator trace, signless trajectory majorant, deletion and source proof structure.
* `ADIABATIC/LINEAR5/STAGE7R/RIGOR.md`, section 3 and Theorem 3.3: unconditional edge SOS, improved trace cap and constants.
* `ADIABATIC/STAGE6-PAPER/main.tex`: screening reference, especially its spatial factorization, edge SOS, marked-source and transfer arguments. The assignment identifies its released article as DOI `10.5281/zenodo.24380128`.
* `ADIABATIC/HUBHOL/RIGOR.md`, sections 1 to 4, and `ADIABATIC/HUBHOL-REVIEW/REVIEW.md`: the on-site Hubbard term and the charge-gap convention.
* `LINEAR5-REVIEW/REVIEW-STAGE5.md`, section 2; `REVIEW-STAGE6.md`; `REVIEW-STAGE7.md`; `REVIEW-STAGE7R.md`; `REVIEW-STAGE7R-DELTA.md`: corrections that are explicitly screened in section 5.

No external paper's formulas or historical priority claims are used. The terms extended Holstein and lattice Frohlich only identify the assignment's density-displacement Hamiltonian. No material, synthesis or measurement claim is made.

## 1. Model, units and the exact atomic reference

### 1.1 Hamiltonian and Gram kernel

**PROVED (conventions and reduction):** put hbar equal to one. On a square torus of side M at least four, let q_i=n_i-1, Delta_i=c_{i,down}c_{i,up}, and take R finite optical branches. The oscillator coordinates obey [x_{jr},p_{ks}]=i delta_{jk}delta_{rs}. Let K_r>0, M_r>0, Omega_r=sqrt(K_r/M_r)>0, real f_r, t>=0, U_c>=0 and real mu. The centered Hamiltonian is

$$H=-t\sum_{\langle ij\rangle,\sigma}(c_{i\sigma}^\dagger c_{j\sigma}+\mathrm{h.c.})+\sum_{j,r}\left(\frac{p_{jr}^2}{2M_r}+\frac{K_rx_{jr}^2}{2}\right)+\sum_{i,j,r}f_r(i-j)q_ix_{jr}+U_c\sum_i(n_{i\uparrow}-\tfrac12)(n_{i\downarrow}-\tfrac12)-\mu\sum_iq_i.$$

The force f_r has units energy/length, K_r energy/length squared, and the Gram matrix has units energy:

$$U_{ij}=\sum_{j',r}\frac{f_r(i-j')f_r(j-j')}{K_r},\qquad U_{ii}=U_{00}=:u.$$

In the uncentered convention H_ep=sum f_r(i-j)n_i x_{jr}, shifting x_{jr} by minus sum_i f_r(i-j)/K_r produces the centered coupling and a linear density term. Translation invariance gives mu_centered=mu_uncentered+sum_j U_{0j}-U_c/2 when the uncentered repulsion is U_c n_up n_down. Additive constants cancel from normalized correlations. There is no restriction on this resulting chemical potential.

**PROVED; CERTIFIED finite folding witness:** for infinite-lattice f_r in l1(Z^2), U(d)=sum_{m,r}f_r(-m)f_r(d-m)/K_r is absolutely summable, with sum_d |U(d)|<=sum_r ||f_r||_1^2/K_r. Periodizing f_r onto the torus periodizes U, not merely its zero-diagonal part:

$$U^{(\Lambda)}(d)=\sum_{c\in M\mathbb Z^2}U(d+c).$$

Absolute convergence justifies rearrangement. If Delta_1=sum_{d!=0}|U(d)|, then Delta_{1,Lambda}<=Delta_1 and U_{00,Lambda}>=U(0)-Delta_1. A dominance condition must specify whether it is finite-volume or infinite-volume. No theorem below is inferred just from dominance.

### 1.2 Lang-Firsov transformation

**PROVED; CERTIFIED finite algebra witness:** write x_{jr}=(a_{jr}+a_{jr}^dagger)/sqrt(2M_r Omega_r) and

$$S=\sum_{i,j,r}\frac{f_r(i-j)}{\sqrt{2K_r\Omega_r}}q_i(a_{jr}^\dagger-a_{jr}).$$

The anti-Hermitian S gives e^S a_{jr} e^{-S}=a_{jr}-sum_i f_r(i-j)q_i/sqrt(2K_r Omega_r). Completing the oscillator square gives the electronic atomic energy, up to a scalar,

$$D_0(q)=-\tfrac12\sum_{ij}U_{ij}q_iq_j+U_c\sum_i(\tfrac14-\tfrac12\mathbf1[n_i=1])-\mu\sum_iq_i.$$

Equivalently, with U_eff=u-U_c,

$$D_0=\frac{U_{\rm eff}}2N_{\rm odd}-\frac12\sum_{i\ne j}U_{ij}q_iq_j-\mu\sum_iq_i+\mathrm{constant}.$$

The hopping c_i^dagger c_j acquires a Weyl displacement with coefficient [f_r(i-m)-f_r(j-m)]/sqrt(2K_r Omega_r) on oscillator (m,r). The on-site Hubbard term is unchanged. It is incorrect to replace the atomic energy by a sum of sites when U_ij!=0 off the diagonal.

### 1.3 Atomic parity gap, including nonlocal interactions

**PROVED; CERTIFIED finite samples:** suppose U_eff>0 and set G=U_eff/2. For arbitrary real off-diagonal U_ij and mu, let E_even be the minimum of D_0 over q_i in {-1,+1} at all sites. Then every q_i in {-1,0,+1} satisfies

$$D_0(q)\ge E_{\rm even}+G\,N_{\rm odd}(q).$$

Proof: replace each zero q_i independently by a symmetric random sign. Every off-diagonal product and the linear chemical-potential term have unchanged expectation. Each replaced site's quadratic term lowers the expected energy by U_eff/2. Therefore some all-even configuration has energy at most D_0(q)-G N_odd. Taking the minimum proves the inequality. In particular the atomic ground energy equals E_even. G is a lower bound on the cost per odd site, not necessarily the exact spectral gap at nonzero mu. In the local model at mu=0 it is attained by making one site odd. The energy of two odd sites above an even configuration is U_eff=2G there.

This argument makes no smallness assumption on Delta_1. It supplies an atomic energy bound but does not supply independent site Gibbs weights or the needed polymer factorization. The gap-to-correlation implication for general nonlocal U remains **OPEN** here.

## 2. The on-site multi-mode theorem

### 2.1 Statement with explicit constants

**Theorem 2.1 (PROVED; constants CERTIFIED).** In section 1.1 set f_r(d)=F_r delta_{d,0}, with arbitrary real F_r, finitely many branches, K_r>0 and Omega_r>0. Define

$$u=\sum_{r=1}^R\frac{F_r^2}{K_r},\quad U_{\rm eff}=u-U_c>0,\quad G=\frac{U_{\rm eff}}2,\quad L=1024,\quad C_U=2^{27},\quad C_G=2^{26}.$$

For every beta>0, every real mu, every square torus of side at least four, and every pair of distinct sites, if beta t^2/U_eff<=2^{-27}, then

$$|\langle\Delta_x^\dagger\Delta_y\rangle_{\rm GC}|\le A q^{d(x,y)},\qquad A=\frac{240}{59},\quad q=\frac5{64}.$$

The finite-volume on-site pair-correlation matrix has absolute row sum at most

$$1+\frac{4Aq}{(1-q)^2}=\frac{512579}{205379}.$$

Under the same on-site pair-order definition and real-source canonical transfer as STAGE7R Theorem 3.3 and the reviewed HUBHOL theorem, at every density rho in (0,2),

$$k_BT_{\rm pair}\le 2^{27}\frac{t^2}{U_{\rm eff}}=2^{26}\frac{t^2}{G}.$$

This excludes long-range or nonexponential order in the on-site electronic pair field in the stated temperature region. It is not a theorem for every possible nonlocal pair wavefunction. Constants are independent of R, every F_r and every positive Omega_r subject to the stated u and U_eff, as well as volume and chemical potential. An infinite number of branches and Omega_r=0 itself are outside the theorem.

### 2.2 Chronological phonon trace and locality

**PROVED:** on a closed Fock trajectory, all oscillator events are time ordered with the latest time on the left, including the two marked pair insertions when present. The normalized on-site oscillator factor is exp(-Phi), where

$$\Phi=\frac14\sum_{x,r}\int_0^\beta\!\int_0^\beta[n_x(s)-n_x(s')]^2K_{r,\beta}(s-s')\,ds\,ds',$$

$$K_{r,\beta}(w)=\frac{U_r\Omega_r}{2}\frac{e^{-\Omega_rw}+e^{-\Omega_r(\beta-w)}}{1-e^{-\beta\Omega_r}},\quad U_r=F_r^2/K_r,\quad 0\le w\le\beta,\quad \int_0^\beta K_{r,\beta}(w)dw=U_r.$$

Apply the corrected infinite-oscillator trace identity branch by branch. Each kernel is nonnegative, so 0<exp(-Phi)<=1; summing positive branch actions introduces no factor R. Disjoint site supports use disjoint oscillator collections. Their atomic electronic weights are also independent. Fermion-even physical bonds and pair sources on disjoint supports commute. Chronological shuffles therefore produce exactly the same signed hard-core spatial gas as the local model, with no imaginary-time factorization within a polymer. The same argument holds for marked pair activities because the closed trajectory includes their occupation jumps. No finite oscillator cutoff is used.

### 2.3 Majorant, gap absorption and cap

**PROVED; CERTIFIED finite CAR/SOS and polynomial witnesses:** strip oscillator factors and take absolute values term by term, rather than assuming a many-body sign gauge. Let B_E be the signless hop matrix on connected physical edge support E, k=|E| and v=|S_E|. As in STAGE6 section 3,

$$|w(E)|,|w_{xy}(E)|\le L^{-k}P_E(L),\qquad P_E(L)=\frac{\operatorname{Tr}_{\rm el}e^{-\beta(D_{0,S_E}-LtB_E)}}{z_a^v}.$$

In the local case D_{0,x}=(U_eff/2)1_odd-mu n_x up to a scalar. The reviewed bond SOS is, for every alpha>0 and g>=0,

$$gB_e\le\alpha(\mathbf1_{\rm odd,u}+\mathbf1_{\rm odd,v})+2g^2\alpha^{-1}I.$$

Its constant is 2, not 1; its proof works in each of the exact bond particle-number sectors N_e in {0,...,4} (three algebraic classes, against the four site cost operators), and its bond gauge is not a global graph gauge. Sum over physical degree at most four, choosing alpha=U_eff/16=G/8. This leaves U_eff/4 per odd site and scalar 32kg^2/U_eff=16kg^2/G. There is no cutoff on g/U_eff or t/U_eff.

Let z=exp(beta mu), w=exp(-beta U_eff/4). Then

$$z_a=1+z^2+2zw^2,\quad z_b=1+z^2+2zw,\quad 5z_a-4z_b=(1-z)^2+2z(2w-1)^2+2zw^2\ge0.$$

Consequently z_b/z_a<=5/4 uniformly in beta and mu. Trace monotonicity of a finite Hermitian matrix exponential follows from ordered eigenvalues, not operator monotonicity of exp. It gives

$$P_E(L)\le(5/4)^v\exp(32k\beta L^2t^2/U_{\rm eff}).$$

The hypothesis beta t^2/U_eff<=2^{-27} yields 32L^2 beta t^2/U_eff<=1/4. Since v<=k+1,

$$|w(E)|,|w_{xy}(E)|\le\frac54\left[\frac{(5/4)e^{1/4}}{L}\right]^k.$$

### 2.4 Spatial count, deletion and canonical source

**PROVED; CERTIFIED constants:** a connected support of k distinct physical bonds containing a fixed root has at most 16^k possibilities, using a deterministic Euler walk of its doubled edges. This count includes cycles. With deletion's factor exp(v), e<3 and e^{1/4}<4/3, the rooted marked sum at distance d is bounded by

$$\frac{15}{4}\sum_{k\ge d}q^k=Aq^d,\qquad q=\frac{16(5/4)\,3(4/3)}{1024}=\frac5{64}.$$

The physical Kotecky-Preiss budget is (15/4)q/(1-q)=75/236<1. The partition function is not replaced by one. The exact deletion ratio for a removed support contributes at most exp(v), which is included in this calculation.

For real-source canonical transfer use the same two site-pair source labels per marked displacement as STAGE6 section 7 and STAGE7R Theorem 3.3. Their operators have norm at most two and preserve local parity. The labelled multigraph has degree at most six, so the rooted count is 36^k. The source disk is |xi|<=tau/(16L)=tau/2^{14}, with tau=1/beta. Its scalar penalty per source label is at most 1/8. The total per-label exponent is therefore still at most 1/4. The source ratio and budget are

$$q_{\rm src}=\frac{36(5/4)\,3(4/3)}{L}=\frac{45}{256},\quad \frac{15}{4}\frac{q_{\rm src}}{1-q_{\rm src}}=\frac{675}{844}<1.$$

All source estimates remain uniform in mu and the marked displacement. The source and unperturbed traces are finite for finitely many modes at positive frequency; free oscillator factors cancel in normalized estimates. Bounded real sources and nearest-neighbor hopping give subextensive boundaries. The reviewed real-source free-energy/Legendre argument therefore applies without a new density hypothesis. It gives the same canonical conclusion, not a claimed finite-volume canonical bound at every fixed particle number. The row sum follows by at most 4d sites at torus distance d>=1 and ||Delta_x||<=1.

## 3. Nonlocal kernels with an exact orthogonal reduction

### 3.1 Common-frequency theorem

**Theorem 3.1 (PROVED; CERTIFIED finite orthogonal example).** Suppose all R branches have a common Omega>0, and on each finite square torus the normalized force rows obey the exact condition

$$\sum_{j,r}\frac{f_r(i-j)f_r(k-j)}{K_r}=u\,\delta_{ik},\qquad u>U_c.$$

Then Theorem 2.1 holds with U_eff=u-U_c and G=U_eff/2, even if individual f_r are spatially nonlocal. The condition is a full Gram identity, not small off-diagonal row sum.

Proof: put Q_{jr}=sqrt(M_r Omega)x_{jr} and P_{jr}=p_{jr}/sqrt(M_r Omega). The oscillator Hamiltonian is (Omega/2)sum(P_{jr}^2+Q_{jr}^2). The vectors h_i(j,r)=f_r(i-j)/sqrt(K_r u) are orthonormal. Complete them to a real orthonormal basis of the R|Lambda| oscillator space and rotate both Q and P by the same basis change. This is a finite canonical orthogonal transformation implemented by a unitary. The oscillator Hamiltonian is invariant. Its coupled part becomes sum_i sqrt(Omega u) q_i Q'_i, with independent bright oscillators and uncoupled dark oscillators. Each bright mode has g_0=sqrt(Omega u/2), hence 2g_0^2/Omega=u. Dark factors cancel from every normalized electronic trace. The electronic pair operators are unchanged by this oscillator-only transformation. Apply Theorem 2.1 with one effective mode per site. This proof is exact at finite volume.

An explicit finite example on M=4 is f=(1/2,1/2,-1/2,1/2) along the horizontal four-site cycle, zero vertically, K=1. Its cyclic autocorrelations are (1,0,0,0). Repeated on each row of the 4 by 4 torus, it has u=1 and obeys the condition. This example is not an assertion of an l1 infinite-lattice kernel with the same four-periodic pattern. For an infinite-lattice l1 kernel, zero off-diagonal infinite autocorrelation implies the finite identity after proper periodization.

### 3.2 Several frequency groups and the necessary limitation

**PROVED:** partition branches into groups of equal frequency. The reduction of section 3.1 can be performed separately in every group if that group's Gram matrix is u_g delta_ik, with u_g>=0. Then there is an on-site bright mode per nonzero group and Theorem 2.1 uses u=sum_g u_g.

It is insufficient that only the total static Gram matrix be diagonal when frequencies differ. The oscillator Hamiltonian is not invariant under a rotation mixing unequal-frequency coordinates, and dynamic cross kernels need not cancel. A finite witness on a four-cycle has K_1=K_2=1, f_1=(1,1/8,0,0), f_2=(1,-1/8,0,0), Omega_1=1 and Omega_2=3. The total static nearest-neighbor Gram entries cancel exactly, but the disjoint-edge dynamic cross exponent at the times in section 4.3 is nonzero. **CERTIFIED:** the outward interval in the replay excludes zero. A theorem for arbitrary such unequal-frequency static cancellations is **OPEN** here, not implied by Theorem 3.1.

## 4. The precise nonlocal obstruction

### 4.1 What fails, and what does not

**PROVED:** the nonlocal Lang-Firsov transformation itself is valid. For a closed trajectory define a_{jr}(s)=sum_i f_r(i-j)q_i(s)/sqrt(K_r). The Gaussian action is (1/4)sum_{j,r}int_0^beta int_0^beta [a_{jr}(s)-a_{jr}(s')]^2 kappa_{r,beta}(s-s') ds ds', where kappa_{r,beta}(w)=(Omega_r/2)[e^{-Omega_r w}+e^{-Omega_r(beta-w)}]/[1-e^{-beta Omega_r}], extended periodically from 0<=w<=beta. This kernel is nonnegative and has integral one. Thus the action is nonnegative and the complete phonon factor remains at most one.

However, disjoint electronic site supports can act on the same oscillator and interact through off-diagonal U_ij. Dropping the globally bounded phonon factor does not restore a product atomic reference. The STAGE6 equality Z=sum_compatible products of connected physical-support activities consequently cannot be imported. Its deletion estimate and the local scalar trace cap have lost a required premise. This is the precise failure before the downstream pair estimate.

### 4.2 Gram-realized electronic counterexample under dominance

**PROVED; CERTIFIED separation of finite traces:** use a four-cycle, U_c=mu=0, U_00=1 and U_01=U_03=1/16, U_02=0. Then Delta_1=1/8<=U_00/4. It is an exact summable Gram kernel: take one branch f_1=(1/4,1/4) with K_1=1 and another f_2=(1,0) with K_2=8/7, along the horizontal direction. The diagonal is 2(1/4)^2+1/(8/7)=1 and each horizontal nearest-neighbor entry is (1/4)^2=1/16. The same local pattern embeds in a square lattice and in a 4 by 4 torus; the electronic reference separates between its horizontal rows.

Consider the exact electronic model D_0=-1/2 sum_i q_i^2-(1/16)sum_i q_i q_{i+1}, beta=1 and hopping amplitude t=1/8. Keep hopping only on e_0={0,1} and e_2={2,3}. Write Z_empty, Z_0, Z_2 and Z_02 for the electronic partition functions with neither edge, each single edge and both edges, respectively. If the STAGE6 disjoint physical-support rule were valid for these two edges, inclusion-exclusion would require

$$Z_{02}Z_{\emptyset}=Z_0Z_2.$$

The replay constructs actual spinful CAR matrices, decomposes by particle numbers per connected component, and evaluates all four traces with rigorous interval matrix exponentials. The difference's interval excludes zero. At t=0 the atomic measure is already nonproduct, since the 1/16 interaction couples these supports. The finite t witness establishes the failed activity identity itself, not just a nonproduct Gibbs measure. It uses at most 16 by 16 blocks and includes all 256 electronic states. Multiplying by other independent square-torus rows cannot repair the failed identity.

This finite trace is an electronic reference check. It is not falsely presented as a finite-phonon computation. The direct finite-frequency witness is separate in section 4.3; together they show that both independent ingredients of the attempted transplant fail.

### 4.3 Finite-frequency oscillator counterexample under dominance

**PROVED; CERTIFIED outward-interval witness:** on a four-cycle take one branch with K=1, f=(1,1/16,0,0), Omega=6/5 and beta=1. This has U_00=257/256, nearest-neighbor U_01=1/16 and Delta_1=1/8<=U_00/4. A hop and return on e_0 occur at s_1=1/5,s_2=3/5; a hop and return on the disjoint e_2 occur at u_1=2/5,u_2=4/5. These are closed legal single-electron trajectories on each edge.

Let g_Omega(w)=cosh[Omega(beta/2-w)]/sinh(beta Omega/2). The coupling-vector inner product of the two oriented edges is U_02-U_03-U_12+U_13=-1/8. With displacement normalization 1/sqrt(2Omega), the Gaussian formula gives

$$\log\frac{P_{02}}{P_0P_2}=\frac{1}{16\Omega}\left[g_\Omega(|s_1-u_1|)-g_\Omega(|s_1-u_2|)-g_\Omega(|s_2-u_1|)+g_\Omega(|s_2-u_2|)\right].$$

Here |s_1-u_1|=|s_2-u_1|=|s_2-u_2|=1/5 and |s_1-u_2|=3/5, so by reflection g(3/5)=g(2/5) the bracket equals g_Omega(1/5)-g_Omega(2/5), which is positive because g decreases on [0,1/2]. Thus P_02>P_0P_2 despite disjoint electronic edges. The replay independently evaluates each oscillator's chronological Gaussian trace, not just this simplified formula, and its difference interval excludes zero. This example is finite frequency and also embeds along square-torus rows.

### 4.4 Why a small fixed spatial row sum does not directly fix the expansion

**PROVED; CERTIFIED scalar witness:** expanding a nonlocal atomic bond as an additional spatial Mayer label produces e^{beta U_ij q_iq_j}-1. In section 4.2, choose q_i=q_j=1 and beta=16; the absolute scalar label is e-1. The replay separates it from e^{-2}/16. A condition on Delta_1/U_00 alone does not control beta Delta_1: beta is unbounded in the regime beta t^2/G small because t may tend to zero. This rules out this proposed beta-independent small-label estimate, not all ways of treating the nonlocal atomic reference.

**OPEN:** constructing a nonproduct reference or resummation that proves an explicit Omega-uniform pair bound for generic summable nonlocal f, including diagonal-dominant examples. No stronger impossibility claim is made.

## 5. Step-by-step screen against reviewed STAGE6/STAGE7R

| Proof ingredient | Local multi-mode / orthogonal reduction | General nonlocal coupling | Review correction respected |
| :-- | :-- | :-- | :-- |
| Lang-Firsov centering and Hubbard term | exact, section 1.2 | exact, section 1.2 | shift mu for uncentered convention; U_eff=u-U_c |
| Infinite oscillator trace | sum of independent actions, section 2.2 | global nonnegative action, section 4.1 | latest-time-left order and corrected BCH sign, STAGE5 review section 2 |
| Disjoint-support spatial gas | independent sites and modes, or rotate exactly | false, sections 4.2 and 4.3 | no unsupported full-graph CAR gauge or imaginary-time factorization |
| Signless majorant | exact after taking termwise absolute values | no independent atomic-support denominator | full three-channel electronic structure; no global positivity assertion |
| Edge gap/SOS | constant 2, arbitrary hopping, section 2.3 | global atomic gap alone is insufficient | STAGE7R corrected SOS, no obsolete propagation cutoff |
| Fugacity cap | exact 5/4 cap, section 2.3 | independent site ratio unavailable | old STAGE6 cap was uniformly bounded, not exponentially unbounded |
| Rooted support count | degree squared: 16^k or 36^k | physical supports do not factor | 88 rooted 3-edge supports require corrected count, not the discarded bound |
| Trace comparison | eigenvalue trace monotonicity | cannot replace nonlocal energy by local trace | no claim that exp is operator monotone |
| Deletion and pair denominator | physical budget 75/236<1 | required hard-core gas not established | divide by exact partition through deletion, not a lower bound of one |
| Canonical transfer | same real-source disk and budgets, section 2.4 | not imported without source locality estimates | no |xi|/t cutoff; no finite-volume ensemble equivalence assertion |
| Frequency limits | bounds uniform for all positive frequencies | no claim about dynamics from static Gram alone | do not identify slow-phonon sum with a bare two-state Hamiltonian |

The STAGE7R second delta review upheld Theorem 3.3 and its canonical transfer but required three ancillary narrative corrections. None of those false narratives is reused here. No Schrieffer-Wolff effective Hamiltonian, pair-weight lower envelope, slow-frequency bare block or universal perturbative no-go is used.

## 6. Certificate scope and resource bound

**CERTIFIED:** `replay.py` uses Fraction arithmetic for every rational matrix/algebraic witness and mpmath.iv at 70 digits for exponentials. Interval matrix exponentials use a Taylor sum with an induced one-norm remainder exp(||A||_1)||A||_1^{N+1}/(N+1)! added to every matrix entry. The trace witnesses require disjoint outward enclosures; no float-converted endpoint decides a certificate inequality. The program checks 20 named guarded mutants, each required to exit 1 with its exact `REJECTED <name>:` prefix. It tests centering, Gram sum, torus folding, the two failed factorizations, the Mayer estimate, orthogonal rows, atomic gap, CAR bond gauge, SOS coefficient, cap and every inherited cluster/source/row-sum constant family.

**NUMERICAL:** a separate four-site local electronic diagonalization at beta=8 and beta t^2/U_eff=2^{-27} reports two Hubbard parameter choices. It is an electronic Omega=infinity check, not a finite-frequency oscillator simulation and not a proof of a thermodynamic statement.

No finite arithmetic certificate is a formal verification of the infinite oscillator identity, graph-general cluster proof or canonical limit; these are analytic arguments with the reviewed dependencies above. The exact interval electronic calculation uses matrices of size at most 16, and the separate floating-point check uses size 256. Mutants run sequentially. The assigned laptop execution requires no cloud, GPU or extensive data. Observed replay exit codes and hashes are recorded in RESULT.md. No other directory is edited.

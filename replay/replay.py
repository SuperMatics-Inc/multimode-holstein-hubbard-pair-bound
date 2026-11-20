#!/usr/bin/env python3
"""EXTHOL replay: local multi-mode theorem and nonlocal factorization obstruction.

Finite witnesses use exact Fraction arithmetic and outward-rounded mpmath.iv
intervals (70 digits). They do not formally verify an infinite-dimensional
trace identity, a cluster expansion, or a thermodynamic limit.

Checks: torus folding, multi-mode Lang-Firsov centering, a finite orthogonal
kernel, nonlocal oscillator and atomic non-factorization, a Mayer bond weight,
the attractive atomic gap, bond CAR/SOS, the 5/4 cap, and all theorem constants.
An electronic 4-site diagonalization is reported separately as NUMERICAL.

Default exits 0; --all-mutants requires every guarded mutant to exit 1 with
'REJECTED <name>:'.
"""

import argparse
import itertools
import subprocess
import sys
from fractions import Fraction as Q
from math import factorial
from pathlib import Path

import mpmath as mp

mp.iv.dps = 70
iv = mp.iv

MUTANTS = (
    "lf-uncentered",
    "ueff-sign",
    "torus-wrap-missing",
    "multimode-sum-wrong",
    "polymer-factor-nonlocal",
    "phonon-factor-nonlocal",
    "mayer-lowtemp",
    "orthogonal-kernel-wrong",
    "charge-gap-double",
    "sos-c1",
    "jw-gauge",
    "cap-54-wrong",
    "scale-too-small",
    "threshold-constant-wrong",
    "physical-degree-wrong",
    "pair-prefactor-wrong",
    "source-disk-wrong",
    "source-degree-wrong",
    "row-shell-wrong",
    "cap-identity-wrong",
)


class Rejected(Exception):
    pass


def require(cond, mutant, msg):
    if not cond:
        raise Rejected(f"{mutant}: {msg}")


# ---------------------------------------------------------------------------
# Exact rational linear algebra & Fock helpers (modeled on HUBHOL/replay.py)
# ---------------------------------------------------------------------------

def is_psd(M):
    """Exact symmetric LDL^T PSD test with Fraction entries."""
    M = [row[:] for row in M]
    n = len(M)
    for k in range(n):
        p = M[k][k]
        if p < 0:
            return False
        if p == 0:
            if any(M[k][j] != 0 for j in range(k + 1, n)):
                return False
            continue
        for i in range(k + 1, n):
            f = M[i][k] / p
            if f:
                for j in range(k + 1, n):
                    M[i][j] -= f * M[k][j]
    return True


def site_n(s, site):
    return ((s >> (2 * site)) & 1) + ((s >> (2 * site + 1)) & 1)


def site_q(s, site):
    return site_n(s, site) - 1


def n_double(s, site):
    return 1 if site_n(s, site) == 2 else 0


COSTS_4 = {
    "n_odd": lambda s, site: 1 if site_n(s, site) == 1 else 0,
    "n_not1": lambda s, site: 0 if site_n(s, site) == 1 else 1,
    "n_occ": lambda s, site: 0 if site_n(s, site) == 0 else 1,
    "n_hole": lambda s, site: 0 if site_n(s, site) == 2 else 1,
}


def gamma_bond(cost, sites=(0, 1)):
    return [sum(cost(s, site) for site in sites) for s in range(1 << (2 * len(sites)))]


def check_bond(mutant):
    """Build 16-state bond hopping T and majorant B; verify JW sign gauge and ||T|| <= 2."""
    T = [[Q(0)] * 16 for _ in range(16)]
    B = [[Q(0)] * 16 for _ in range(16)]
    for s in range(16):
        for spin in (0, 1):
            for src, dst in ((spin, 2 + spin), (2 + spin, spin)):
                if ((s >> src) & 1) and not ((s >> dst) & 1):
                    lo, hi = (src, dst) if src < dst else (dst, src)
                    between = sum((s >> k) & 1 for k in range(lo + 1, hi))
                    sgn = -1 if (between % 2) else 1
                    s2 = s ^ (1 << src) ^ (1 << dst)
                    T[s2][s] += sgn
                    B[s2][s] += 1
    require(all(T[i][j] == T[j][i] for i in range(16) for j in range(16)), mutant, "T not symmetric")
    require(all(abs(T[i][j]) == B[i][j] for i in range(16) for j in range(16)), mutant, "|T| != B")
    require(all(B[i][j] == B[j][i] for i in range(16) for j in range(16)), mutant, "B not symmetric")

    # The one-up, one-down sector is K_{2,2}.
    # Build T_e^2 and verify 4*I - T^2 is PSD (hence ||T|| <= 2)
    T2 = [[sum(T[i][k] * T[k][j] for k in range(16)) for j in range(16)] for i in range(16)]
    M4 = [[(Q(4) if i == j else Q(0)) - T2[i][j] for j in range(16)] for i in range(16)]
    require(is_psd(M4), mutant, "4I - T_e^2 not PSD")

    # Bits: up0=0, down0=1, up1=2, down1=3.
    # The legal cycle 3 -> 9 -> 12 -> 6 -> 3 has two negative CAR hops.
    # Flipping exactly one hop sign frustrates it.
    k22_states = (3, 9, 12, 6)
    gauge_prod = 1
    for idx_step, (a, b) in enumerate(zip(k22_states, k22_states[1:] + k22_states[:1])):
        val = T[b][a]
        require(abs(val) == 1, mutant, "K22 step missing")
        step_sgn = -int(val) if (mutant == "jw-gauge" and idx_step == 0) else int(val)
        gauge_prod *= step_sgn
    # The cycle product is +1, allowing a diagonal sign gauge on this component.
    require(gauge_prod == 1, mutant, f"K22 cycle gauge product {gauge_prod} != 1")
    return T, B


def check_edge_sos(B, mutant):
    """Certify Gamma + c_sos * x^2 I - x B_e >= 0 for c_sos = 2 (and reject c_sos = 1)."""
    c_sos = Q(1) if mutant == "sos-c1" else Q(2)
    x_samples = [Q(-2), Q(-1), Q(-1, 2), Q(0), Q(1, 2), Q(1), Q(2)]
    for name, cost in COSTS_4.items():
        G_diag = gamma_bond(cost)
        for x in x_samples:
            M = [[-x * B[i][j] for j in range(16)] for i in range(16)]
            for i in range(16):
                M[i][i] += G_diag[i] + c_sos * x * x
            require(is_psd(M), mutant, f"Edge SOS failed for {name} at x={x}, c_sos={c_sos}")


# ---------------------------------------------------------------------------
# Interval matrix exponential via exact Taylor series + rigorous remainder
# ---------------------------------------------------------------------------

def iv_mat_exp_neg(H_iv, beta_iv, terms=26):
    """Compute rigorous enclosure of exp(-beta * H) via Taylor sum + remainder bound."""
    n = len(H_iv)
    # Compute A = -beta * H and bound ||A||_1
    A = [[(-beta_iv) * H_iv[i][j] for j in range(n)] for i in range(n)]
    norm1_hi = iv.mpf(0)
    for j in range(n):
        col_sum = iv.mpf(0)
        for i in range(n):
            col_sum += abs(A[i][j])
        if col_sum.b > norm1_hi:
            norm1_hi = col_sum.b
    norm1_iv = iv.mpf(norm1_hi)

    # S = sum_{k=0}^{terms-1} A^k / k!
    S = [[iv.mpf(1) if i == j else iv.mpf(0) for j in range(n)] for i in range(n)]
    P = [[iv.mpf(1) if i == j else iv.mpf(0) for j in range(n)] for i in range(n)]
    for k in range(1, terms):
        # P = P @ A / k
        P_next = [[iv.mpf(0)] * n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                s = iv.mpf(0)
                for l in range(n):
                    if P[i][l] != 0 and A[l][j] != 0:
                        s += P[i][l] * A[l][j]
                P_next[i][j] = s / iv.mpf(k)
        P = P_next
        for i in range(n):
            for j in range(n):
                S[i][j] += P[i][j]

    # Remainder bound: ||A||_1^{terms} / terms! * exp(||A||_1)
    rem_bound = (norm1_iv ** terms) / iv.mpf(factorial(terms)) * mp.iv.exp(norm1_iv)
    err = iv.mpf([-1, 1]) * rem_bound.b
    for i in range(n):
        for j in range(n):
            S[i][j] += err
    return S


def iv_trace_exp_neg(H_q, beta_q, terms=26):
    """Exact-rational matrix H_q -> rigorous interval Tr(exp(-beta_q * H_q))."""
    n = len(H_q)
    H_iv = [[iv.mpf(H_q[i][j].numerator) / iv.mpf(H_q[i][j].denominator) for j in range(n)] for i in range(n)]
    beta_iv = iv.mpf(beta_q.numerator) / iv.mpf(beta_q.denominator)
    E = iv_mat_exp_neg(H_iv, beta_iv, terms=terms)
    tr = iv.mpf(0)
    for i in range(n):
        tr += E[i][i]
    return tr


# ---------------------------------------------------------------------------
# Main certification suite
# ---------------------------------------------------------------------------

def check_torus_folding(mutant):
    """Certify Lemma 1.1: U_{0j}^{(Lambda)} = sum_c U_{0, j+c}^{(inf)} and Delta_{1,Lambda} <= Delta_{1,inf}."""
    # Exact rational 1D/2D model with finite-support f on Z: f(0)=1, f(1)=1/8, f(-1)=-1/16, f(3)=1/32
    f_inf = {0: Q(1), 1: Q(1, 8), -1: Q(-1, 16), 3: Q(1, 32)}
    K = Q(1)
    # Compute U_{0, d}^{(inf)} = sum_m f_inf(-m) * f_inf(d - m)
    U_inf = {}
    for a, fa in f_inf.items():      # a = -m
        for b, fb in f_inf.items():  # b = d - m = d + a  => d = b - a
            d = b - a
            U_inf[d] = U_inf.get(d, Q(0)) + fa * fb / K

    M = 4
    # Fold f onto Z/MZ: f_M(r) = sum_{a = r mod M} f_inf(a)
    f_M = {r: Q(0) for r in range(M)}
    for k, fk in f_inf.items():
        f_M[k % M] += fk

    # Direct torus sum: U_{0,j}^{(M)} = sum_{m=0}^{M-1} f_M((-m)%M) * f_M((j-m)%M) / K
    U_M_direct = {j: Q(0) for j in range(M)}
    for j in range(M):
        for m in range(M):
            U_M_direct[j] += f_M[(-m) % M] * f_M[(j - m) % M] / K

    # Folded sum from U_inf: sum_{c in M*Z} U_inf[j + c]
    U_M_folded = {j: Q(0) for j in range(M)}
    for d, u in U_inf.items():
        if mutant == "torus-wrap-missing" and d != d % M:
            continue
        U_M_folded[d % M] += u

    for j in range(M):
        require(U_M_direct[j] == U_M_folded[j], mutant, f"Torus folding mismatch at j={j}")

    Delta1_inf = sum(abs(u) for d, u in U_inf.items() if d != 0)
    Delta1_M = sum(abs(U_M_direct[j]) for j in range(1, M))
    require(Delta1_M <= Delta1_inf, mutant, "Delta_{1,Lambda} > Delta_{1,inf}")
    require(U_M_direct[0] >= U_inf[0] - Delta1_inf, mutant, "U00^{(Lambda)} bound failed")


def check_lf_and_multimode(mutant):
    """Certify Lemma 1.2: exact multi-mode nonlocal Lang-Firsov centering & U_eff / mu shift."""
    # Two modes r=1,2 on a 2-site ring
    K1, K2 = Q(3, 2), Q(5, 3)
    f1 = [Q(1), Q(1, 4)]
    f2 = [Q(2, 5), Q(-1, 5)]
    Uc = Q(7, 5)
    mu_raw = Q(3, 11)

    if mutant == "multimode-sum-wrong":
        # Wrongly take max_r instead of sum_r
        U00 = max(sum(x * x for x in f1) / K1, sum(x * x for x in f2) / K2)
        U01 = max(2 * f1[0] * f1[1] / K1, 2 * f2[0] * f2[1] / K2)
    else:
        U00 = sum(x * x for x in f1) / K1 + sum(x * x for x in f2) / K2
        U01 = 2 * f1[0] * f1[1] / K1 + 2 * f2[0] * f2[1] / K2

    Ueff = (Uc - U00) if mutant == "ueff-sign" else (U00 - Uc)
    if mutant == "lf-uncentered":
        mu_c = mu_raw - Uc / 2
    else:
        mu_c = mu_raw - Uc / 2 + U00 + U01

    diffs = set()
    for n0 in (0, 1, 2):
        for n1 in (0, 1, 2):
            p0 = 1 if n0 == 2 else 0
            p1 = 1 if n1 == 2 else 0
            # Note: compute true sum_r U_r for D_raw check
            U00_true = sum(x * x for x in f1) / K1 + sum(x * x for x in f2) / K2
            U01_true = 2 * f1[0] * f1[1] / K1 + 2 * f2[0] * f2[1] / K2
            D_raw = (
                Uc * (p0 + p1)
                - mu_raw * (n0 + n1)
                - Q(1, 2) * (U00_true * n0 * n0 + 2 * U01_true * n0 * n1 + U00_true * n1 * n1)
            )
            q0, q1 = n0 - 1, n1 - 1
            D_cen = -Ueff * Q(1, 2) * (q0 * q0 + q1 * q1) - U01 * q0 * q1 - mu_c * (q0 + q1)
            diffs.add(D_raw - D_cen)
    require(len(diffs) == 1, mutant, f"LF centering non-constant diff: {diffs}")


def check_barker4_and_phonon_counterexample(mutant):
    """Certify Barker-4 orthogonal equivalence & finite-Omega phonon non-factorization."""
    # 1. Barker-4 on 4-site ring: f = (1/2, 1/2, -1/2, 1/2) -> orthogonal F^T F = I
    f_b = [Q(1, 2), Q(1, 2), Q(-1, 2), Q(1, 2)]
    if mutant == "orthogonal-kernel-wrong":
        f_b[2] = Q(1, 2)
    for j in range(4):
        u0j = sum(f_b[(-m) % 4] * f_b[(j - m) % 4] for m in range(4))
        require(u0j == (Q(1) if j == 0 else Q(0)), mutant, f"Barker-4 failed at j={j}")

    # Finite-Omega non-factorization with f=(1,1/16,0,0), K=1, Omega=6/5, beta=1.
    # Hop e0={0,1} at s1=1/5, s2=3/5; hop e2={2,3} at u1=2/5, u2=4/5
    # Exact Gaussian formula: exp(-0.5 * sum_{k,l} c_k c_l * cosh(Omega*(beta/2 - |t_k - t_l|)) / sinh(beta*Omega/2))
    Omega = iv.mpf(6) / iv.mpf(5)
    beta = iv.mpf(1)
    s1, u1, s2, u2 = iv.mpf(1) / 5, iv.mpf(2) / 5, iv.mpf(3) / 5, iv.mpf(4) / 5

    def G_ph(w):
        a = Omega * (beta / 2 - w)
        b = beta * Omega / 2
        return (mp.iv.exp(a) + mp.iv.exp(-a)) / (mp.iv.exp(b) - mp.iv.exp(-b))

    def ph_factor(events):
        phi = iv.mpf(0)
        for tk, ck in events:
            for tl, cl in events:
                phi += iv.mpf("0.5") * ck * cl * G_ph(abs(tk - tl))
        return mp.iv.exp(-phi)

    f_nn = [iv.mpf(1), iv.mpf(1) / 16, iv.mpf(0), iv.mpf(0)]
    P_both = iv.mpf(1)
    P_e0 = iv.mpf(1)
    P_e2 = iv.mpf(1)
    for m in range(4):
        # g/Omega = sqrt(1/(2*Omega)) * f  => c = (g(u-m)-g(v-m))/Omega => c^2 = (f_u - f_v)^2 / (2*Omega)
        pref = iv.mpf(1) / mp.iv.sqrt(iv.mpf(2) * Omega)
        c0 = pref * (f_nn[(0 - m) % 4] - f_nn[(1 - m) % 4])
        c2 = pref * (f_nn[(2 - m) % 4] - f_nn[(3 - m) % 4])
        P_both *= ph_factor([(s1, c0), (u1, c2), (s2, -c0), (u2, -c2)])
        P_e0 *= ph_factor([(s1, c0), (s2, -c0)])
        P_e2 *= ph_factor([(u1, c2), (u2, -c2)])
    require(Q(1, 8) <= Q(257, 256) / 4, mutant, "Dominance condition failed")

    diff_iv = P_both - P_e0 * P_e2
    contains_zero = (diff_iv.a <= 0) and (0 <= diff_iv.b)
    if mutant == "phonon-factor-nonlocal":
        # Falsely assert factorization holds (0 in diff_iv)
        require(contains_zero, mutant, "Phonon trace failed to factorize across disjoint bonds")
    else:
        require(not contains_zero, mutant, "Failed to separate P_both from P_e0*P_e2")

    # 3. Also certify 2-mode opposite-U01 counterexample (total static U01=0, Omega1=1, Omega2=3)
    Omega1, Omega2 = iv.mpf(1), iv.mpf(3)
    # Mode 1: f1=[1,1/8,0,0], K1=1; mode 2: f2=[1,-1/8,0,0], K2=1.
    # Each nearest-neighbor U01 is +/-1/8. Their static sum is zero.
    # The disjoint-edge inner product is 2*U02-U03-U01 = -2*U01.
    def cross_kernel(Om):
        def g(w):
            a = Om * (beta / 2 - w)
            b = beta * Om / 2
            return (mp.iv.exp(a) + mp.iv.exp(-a)) / (mp.iv.exp(b) - mp.iv.exp(-b))
        return (g(abs(s1 - u1)) - g(abs(s1 - u2)) - g(abs(s2 - u1)) + g(abs(s2 - u2))) / (iv.mpf(2) * Om)

    # Dynamic cross exponent magnitude is (1/4)*(cross_kernel(O1)-cross_kernel(O2)).
    dyn_cross = (iv.mpf(1) / 4) * (cross_kernel(Omega1) - cross_kernel(Omega2))
    require(not (dyn_cross.a <= 0 <= dyn_cross.b), mutant, "2-mode canceling static U01 dynamic cross failed")
    if mutant is None:
        print(f"CERTIFIED: oscillator P_both-P_e0*P_e2 = {diff_iv}")
        print(f"CERTIFIED: unequal-frequency static-cancellation cross exponent = {dyn_cross}")


def check_atomic_nonfactorization_and_mayer(mutant):
    """Certify Part (b): w({e0,e2}) != w({e0})w({e2}) at Omega=inf even with Delta_1 = U00/8 <= U00/4,
    and Mayer bond weight breakdown at low T."""
    # 4-site ring: U00 = 1, U01 = 1/16 (so Delta_1 = 1/8 <= U00/4), Uc = 0, mu = 0, beta = 1, t = 1/8.
    # For edge subset E in {(), (e0,), (e2,), (e0,e2)}:
    # Z(E) = Tr(exp(-beta * H_E)) on 256 states -> split by conserved quantities!
    # Conserved quantities for E={(0,1)}: n_0+n_1 (up/dn separate? No, D0 preserves every n_{i,sig} only when t=0;
    # for E={(0,1)}, site occupations (n_{2,up},n_{2,dn},n_{3,up},n_{3,dn}) and (n_{0,up}+n_{1,up}, n_{0,dn}+n_{1,dn}) are conserved!
    # Max block size for E={(0,1)} is 4x4; for E={(0,1),(2,3)} is 16x16 (still tiny and exact!).
    U00 = Q(1)
    U01 = Q(1, 16)
    beta = Q(1)
    t = Q(1, 8)
    # Exact Gram realization: f1=(1/4,1/4), K1=1; f2=(1,0), K2=8/7.
    require(2 * Q(1, 4) ** 2 + Q(1) / Q(8, 7) == U00,
            mutant, "Gram diagonal mismatch")
    require(Q(1, 4) ** 2 == U01, mutant, "Gram off-diagonal mismatch")
    require(2 * U01 == U00 / 8 and 2 * U01 <= U00 / 4,
            mutant, "Atomic witness dominance mismatch")

    def block_partition(edges):
        # Partition 256 states by connected components of edges + spin counts on each component
        parent = list(range(4))
        for u, v in edges:
            # merge u, v
            pu, pv = u, v
            while parent[pu] != pu:
                pu = parent[pu]
            while parent[pv] != pv:
                pv = parent[pv]
            if pu != pv:
                parent[pv] = pu
        comps = {}
        for i in range(4):
            p = i
            while parent[p] != p:
                p = parent[p]
            comps.setdefault(p, []).append(i)
        comp_list = list(comps.values())

        def key_of(s):
            return tuple(
                (sum((s >> (2 * i)) & 1 for i in c), sum((s >> (2 * i + 1)) & 1 for i in c))
                for c in comp_list
            )

        buckets = {}
        for s in range(256):
            buckets.setdefault(key_of(s), []).append(s)

        Z_tot = iv.mpf(0)
        for states in buckets.values():
            m = len(states)
            idx = {s: k for k, s in enumerate(states)}
            H = [[Q(0)] * m for _ in range(m)]
            for k, s in enumerate(states):
                qs = [site_q(s, i) for i in range(4)]
                d0 = -Q(1, 2) * U00 * sum(q * q for q in qs) - U01 * sum(qs[i] * qs[(i + 1) % 4] for i in range(4))
                H[k][k] = d0
                for u, v in edges:
                    for sig in (0, 1):
                        bu, bv = 2 * u + sig, 2 * v + sig
                        for src, dst in ((bu, bv), (bv, bu)):
                            if ((s >> src) & 1) and not ((s >> dst) & 1):
                                lo, hi = (src, dst) if src < dst else (dst, src)
                                sgn = -1 if (sum((s >> b) & 1 for b in range(lo + 1, hi)) % 2) else 1
                                s2 = s ^ (1 << src) ^ (1 << dst)
                                H[idx[s2]][k] -= t * sgn
            Z_tot += iv_trace_exp_neg(H, beta, terms=22)
        return Z_tot

    Z_empty = block_partition(())
    Z_e0 = block_partition(((0, 1),))
    Z_e2 = block_partition(((2, 3),))
    Z_both = block_partition(((0, 1), (2, 3)))

    # Polymer inclusion-exclusion activities (normalized by Z_empty):
    # w({e0}) = Z_e0/Z_empty - 1; w({e2}) = Z_e2/Z_empty - 1;
    # w({e0,e2}) = Z_both/Z_empty - Z_e0/Z_empty - Z_e2/Z_empty + 1 = w({e0})*w({e2}) iff Z_both*Z_empty == Z_e0*Z_e2!
    diff = Z_both * Z_empty - Z_e0 * Z_e2
    if mutant == "polymer-factor-nonlocal":
        require(iv.mpf(0) in diff, mutant, "Atomic polymer weight failed to factorize across disjoint edges")
    else:
        require(not (iv.mpf(0) in diff), mutant, "Failed to separate Z_both*Z_empty from Z_e0*Z_e2")
    if mutant is None:
        print(f"CERTIFIED: atomic Z_both*Z_empty-Z_e0*Z_e2 = {diff}")

    # Mayer low-temperature failure: for U01 = 1/16, beta * U01 = 1 -> e^{beta*U01} - 1 = e - 1 > e^{-2}/16
    beta_low = iv.mpf(16)
    mayer_wt = mp.iv.exp(beta_low * iv.mpf(1) / 16) - iv.mpf(1)
    kp_thresh = mp.iv.exp(iv.mpf(-2)) / iv.mpf(16)
    if mutant == "mayer-lowtemp":
        require(mayer_wt.b <= kp_thresh.a, mutant, "Mayer bond weight exceeds KP threshold")
    else:
        require(mayer_wt.a > kp_thresh.b, mutant, "Mayer threshold check failed")


def check_atomic_gaps(mutant):
    """Check the attractive odd-site gap against exact 4-site spectra."""
    # 1. Attractive U_eff > 0: test 3^4 = 81 configurations for arbitrary U_ij (including large |U_01| > U00) and mu
    for U00, Uc, U01, U02, mu in [
        (Q(1), Q(0), Q(0), Q(0), Q(0)),
        (Q(1), Q(0), Q(1, 4), Q(0), Q(0)),
        (Q(1), Q(0), Q(-3, 2), Q(1, 3), Q(2, 5)),
        (Q(3, 2), Q(1, 2), Q(-2), Q(-1, 2), Q(7, 4)),
    ]:
        Ueff = U00 - Uc
        E_even = []
        E_odd = []
        all_samples = []
        for qs in itertools.product((-1, 0, 1), repeat=4):
            E = (
                -Ueff * Q(1, 2) * sum(q * q for q in qs)
                - U01 * sum(qs[i] * qs[(i + 1) % 4] for i in range(4))
                - U02 * (qs[0] * qs[2] + qs[1] * qs[3])
                - mu * sum(qs)
            )
            all_samples.append((E, qs.count(0)))
            if 0 in qs:
                E_odd.append(E)
            else:
                E_even.append(E)
        gap_odd = min(E_odd) - min(E_even)
        claimed_gap = Ueff if mutant == "charge-gap-double" else Ueff / 2
        require(gap_odd >= claimed_gap, mutant,
                f"Odd-sector gap {gap_odd} < claimed={claimed_gap}")
        require(all(E >= min(E_even) + claimed_gap * n for E, n in all_samples),
                mutant, "Per-odd-site gap inequality failed")

    # Witness attaining equality: U00=1, Uc=0, U01=-1/2, U02=0, mu=1 -> odd gap == 1/2
    Ueff = Q(1)
    U01 = Q(-1, 2)
    mu = Q(1)
    E_even = []
    E_odd = []
    for qs in itertools.product((-1, 0, 1), repeat=4):
        E = -Ueff * Q(1, 2) * sum(q * q for q in qs) - U01 * sum(qs[i] * qs[(i + 1) % 4] for i in range(4)) - mu * sum(qs)
        (E_odd if 0 in qs else E_even).append(E)
    require(min(E_odd) - min(E_even) == Q(1, 2), mutant, "Equality witness for G_A failed")



def check_caps_and_constants(mutant):
    """Exact cap identity and every finite constant in Theorem 2.1."""
    # 5*z_a-4*z_b equals (1-z)^2+2*z*(2*w-1)^2+2*z*w^2.
    last = 4 if mutant == "cap-identity-wrong" else 2
    for z in (Q(0), Q(1), Q(2)):
        for w in (Q(0), Q(1), Q(2)):
            za = 1 + 2 * z * w ** 2 + z ** 2
            zb = 1 + 2 * z * w + z ** 2
            rhs = (1 - z) ** 2 + 2 * z * (2 * w - 1) ** 2 + last * z * w ** 2
            require(5 * za - 4 * zb == rhs, mutant, "Cap polynomial identity mismatch")
    c54 = Q(1) if mutant == "cap-54-wrong" else Q(5, 4)
    for num in range(100):
        v = iv.mpf([num, num + 1]) / 100
        polynomial = iv.mpf(c54.numerator) / c54.denominator * (1 + v ** 2) - (1 + v)
        require(polynomial.a >= 0, mutant,
                f"5/4 cap failed on [{num}/100,{(num+1)}/100]")

    e_upper = sum((Q(1, factorial(n)) for n in range(6)), Q(0)) + Q(7, 2160)
    require(e_upper == Q(1175, 432) and e_upper < Q(49, 18) < 3,
            mutant, "Exponential rational enclosure mismatch")
    require(Q(4, 3) ** 4 > 3, mutant, "Fourth-root enclosure mismatch")
    require(mp.iv.exp(iv.mpf(1) / 4).b < (iv.mpf(4) / 3).a,
            mutant, "exp(1/4) enclosure failed")

    L = 2 ** 10
    C = Q(2 ** 26) if mutant == "threshold-constant-wrong" else Q(2 ** 27)
    epsilon = Q(1, 2 ** 27)
    if mutant == "scale-too-small":
        epsilon *= 2
    require(C == 128 * L ** 2, mutant, "Energy threshold constant mismatch")
    require(C / 2 == Q(67108864), mutant, "Charge-gap threshold conversion mismatch")
    delta = 32 * L ** 2 * epsilon
    require(delta <= Q(1, 4), mutant, f"Physical exponent {delta} exceeds 1/4")
    count = Q(8) if mutant == "physical-degree-wrong" else Q(16)
    prefactor = Q(5, 4) if mutant == "pair-prefactor-wrong" else Q(15, 4)
    q = count * Q(5, 4) * Q(3) * Q(4, 3) / L
    A = prefactor / (1 - q)
    budget = prefactor * q / (1 - q)
    require(q == Q(5, 64) and A == Q(240, 59) and budget == Q(75, 236) < 1,
            mutant, "Physical cluster constants mismatch")
    shell = 2 if mutant == "row-shell-wrong" else 4
    row_sum = 1 + shell * A * q / (1 - q) ** 2
    require(row_sum == Q(512579, 205379), mutant, "Pair row sum mismatch")
    source_prefactor = Q(1, (8 if mutant == "source-disk-wrong" else 16) * L)
    source_penalty = 2 * L * source_prefactor
    source_count = Q(16) if mutant == "source-degree-wrong" else Q(36)
    source_q = source_count * Q(5, 4) * 3 * Q(4, 3) / L
    source_budget = Q(15, 4) * source_q / (1 - source_q)
    require(source_prefactor == Q(1, 2 ** 14) and source_penalty == Q(1, 8),
            mutant, "Source disk constants mismatch")
    require(source_q == Q(45, 256) and source_budget == Q(675, 844) < 1,
            mutant, "Source cluster constants mismatch")


def numerical_crosschecks(verbose=False):
    """NUMERICAL cross-checks: 4-site electronic Omega=infinity limit, not proofs."""
    import numpy as np

    N = 4
    dim = 4 ** N
    n_op = np.zeros((N, dim))
    for s in range(dim):
        for i in range(N):
            n_op[i, s] = ((s >> (2 * i)) & 1) + ((s >> (2 * i + 1)) & 1)
    q_op = n_op - 1.0

    def el_hop(u, v):
        M = np.zeros((dim, dim))
        for s in range(dim):
            for sig in (0, 1):
                bu, bv = 2 * u + sig, 2 * v + sig
                if ((s >> bv) & 1) and not ((s >> bu) & 1):
                    lo, hi = min(bu, bv), max(bu, bv)
                    sgn = -1.0 if (sum((s >> k) & 1 for k in range(lo + 1, hi)) % 2) else 1.0
                    s2 = s ^ (1 << bv) ^ (1 << bu)
                    M[s2, s] += -sgn
                    M[s, s2] += -sgn
        return M

    T_sum = el_hop(0, 1) + el_hop(1, 2) + el_hop(2, 3) + el_hop(3, 0)
    P02 = np.zeros((dim, dim))
    for s in range(dim):
        if ((s >> 4) & 1) and ((s >> 5) & 1) and not ((s >> 0) & 1) and not ((s >> 1) & 1):
            P02[s ^ (1 << 4) ^ (1 << 5) ^ (1 << 0) ^ (1 << 1), s] = 1.0

    max_ratio = 0.0
    for U00, Uc in [(1.0, 0.0), (1.0, 0.4)]:
        U01, mu = 0.0, 0.2
        A, q = 240.0 / 59.0, 5.0 / 64.0
        Ueff = U00 - Uc
        U_mat = np.zeros((N, N))
        for i in range(N):
            U_mat[i, i] = Ueff
            U_mat[i, (i + 1) % N] = U01
            U_mat[(i + 1) % N, i] = U01
        D0 = np.zeros(dim)
        for i in range(N):
            for j in range(N):
                D0 -= 0.5 * U_mat[i, j] * q_op[i] * q_op[j]
        D0 -= mu * np.sum(q_op, axis=0)
        beta = 8.0
        t = np.sqrt(Ueff * 2 ** (-27) / beta)
        H = np.diag(D0) + t * T_sum
        ev, U = np.linalg.eigh(H)
        ev -= ev[0]
        w = np.exp(-beta * ev)
        corr = abs(np.sum(w * np.diag(U.T @ P02 @ U)) / np.sum(w))
        bnd = A * (q ** 2)
        if corr > bnd:
            raise AssertionError(f"NUMERICAL 4-site check exceeded bound: {corr} > {bnd}")
        max_ratio = max(max_ratio, corr / bnd)

    if verbose:
        print(f"NUMERICAL: 4-site local electronic Omega=inf diagonalization: max corr/bound = {max_ratio:.3e} <= 1")


def certificate(mutant=None, verbose=False):
    _, B = check_bond(mutant)
    check_edge_sos(B, mutant)
    check_torus_folding(mutant)
    check_lf_and_multimode(mutant)
    check_barker4_and_phonon_counterexample(mutant)
    check_atomic_nonfactorization_and_mayer(mutant)
    check_atomic_gaps(mutant)
    check_caps_and_constants(mutant)
    if mutant is None:
        numerical_crosschecks(verbose=verbose)
        if verbose:
            print("CERTIFIED: all EXTHOL exact rational and outward-interval checks passed")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mutant", choices=MUTANTS, help="Run a single guarded mutant (must exit 1)")
    parser.add_argument("--all-mutants", action="store_true", help="Verify all guarded mutants exit 1")
    args = parser.parse_args()

    if args.mutant:
        try:
            certificate(mutant=args.mutant, verbose=False)
        except Rejected as exc:
            print(f"REJECTED {args.mutant}: {exc}")
            return 1
        print(f"UNCAUGHT {args.mutant}")
        return 0

    try:
        certificate(verbose=True)
    except Rejected as exc:
        print(f"BASELINE FAILURE: {exc}")
        return 1

    if not args.all_mutants:
        return 0

    script = str(Path(__file__).resolve())
    for mutant in MUTANTS:
        child = subprocess.run(
            [sys.executable, "-B", script, "--mutant", mutant],
            capture_output=True,
            text=True,
            check=False,
        )
        print(child.stdout.strip())
        if child.stderr:
            print(child.stderr.strip())
        print(f"mutant {mutant}: exit {child.returncode}")
        if child.returncode != 1 or not child.stdout.startswith(f"REJECTED {mutant}:"):
            return 1
    print(f"CERTIFIED: all {len(MUTANTS)} mutants rejected with exit 1")
    return 0


if __name__ == "__main__":
    sys.exit(main())

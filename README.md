# Numerical Simulation Pipeline — Coupled RLC Circuit

> From-scratch implementation of ODE solvers, interpolation methods,
> numerical integration, and nonlinear root-finding in MATLAB,
> applied to a coupled inductor circuit with nonlinear mutual inductance.

---

## Overview

This project implements a complete **numerical simulation pipeline**
for a third-order ODE system. Every algorithm is built without
high-level solver abstractions (no `ode45`, no `interp1`) to demonstrate
deep understanding of the numerical foundations underlying
scientific computing and data engineering.

| Part | Topic | Core concept |
|------|-------|-------------|
| 1 | ODE integration | Time-series state simulation, solver stability |
| 2 | Interpolation & approximation | Sparse data reconstruction, feature engineering |
| 3 | Numerical integration | Signal aggregation, energy metering |
| 4 | Nonlinear root-finding | Parameter optimization, model calibration |

---

## System Model

A coupled inductor circuit described by a 3-variable state-space ODE
system derived from Kirchhoff's laws:
di₁/dt = [-R₁/M · i₁ + R₂/L₂ · i₂ - 1/M · uC + 1/M · e] / D₁
di₂/dt = [-R₁/L₁ · i₁ + R₂/M  · i₂ - 1/L₁ · uC + 1/L₁ · e] / D₂
duC/dt = i₁ / C
where D₁ = L₁/M - M/L₂,  D₂ = M/L₁ - L₂/M

**Circuit parameters:**

| Symbol | Value | Description |
|--------|-------|-------------|
| R₁ | 0.1 Ω | Primary resistance |
| R₂ | 10 Ω | Secondary resistance |
| C | 0.5 F | Capacitance |
| L₁ | 3 H | Primary inductance |
| L₂ | 5 H | Secondary inductance |
| M | 0.8 H | Mutual inductance (linear) / M(u₁) nonlinear |

**Forcing signals:**

```matlab
e1 = @(t) (mod(t,3) < 1.5) * 120;          % square wave, T=3s
e2 = @(t) 240 * sin(t);                     % low-frequency sine
e3 = @(t) 210 * sin(2*pi*5*t);             % 5 Hz sine
e4 = @(t) 120 * sin(2*pi*50*t);            % 50 Hz sine
```

---

## Repository Structure
.
├── czesc1.m        # Part 1 — ODE solvers (Euler & Improved Euler)
├── czesc2.m        # Part 2 — Nonlinear inductance interpolation
├── czesc3.m        # Part 3 — Energy & average power integration
├── czesc4.m        # Part 4 — Frequency root-finding
├── Raport.pdf      # Full experiment report with plots and analysis
└── README.md

---

## Part 1 — ODE Solvers

Two explicit integrators implemented from scratch and benchmarked
across four forcing signals with deliberately varied step sizes
to expose stability boundaries.

**Forward Euler:**
y_{n+1} = y_n + h · f(t_n, y_n)

**Improved Euler (Heun / Midpoint method):**
k₁ = f(t_n, y_n)
y* = y_n + (h/2) · k₁
y_{n+1} = y_n + h · f(t_n + h/2, y*)

**Stability comparison:**

| Forcing signal | Step h | Forward Euler | Improved Euler |
|----------------|--------|--------------|----------------|
| Square wave 3s | 0.1 s | Overestimates peaks | Stable ✓ |
| 240·sin(t) | 0.5 s | Diverges over time | Stable ✓ |
| 210·sin(2π·5t) | 0.1 s | Numerical underflow (10⁻¹²) | Stable ✓ |
| 120·sin(2π·50t) | 0.009 s | Unstable oscillations | Stable ✓ |

Improved Euler achieves **second-order accuracy** at the same step size
as Forward Euler — analogous to why higher-order methods matter in
production time-series pipelines when computational cost per step is fixed.

---

## Part 2 — Nonlinear Inductance: Sparse Data Reconstruction

Mutual inductance M is known only at 8 measured operating points
(20–300 V range). Four methods reconstruct M(u₁) continuously:

| Method | Type | Data points |
|--------|------|-------------|
| Lagrange polynomial | Interpolation | 8 non-uniform |
| Cubic B-spline | Interpolation | 10 uniform (manually resampled) |
| Least-squares deg. 3 | Approximation | 8 non-uniform |
| Least-squares deg. 5 | Approximation | 8 non-uniform (scaled) |

**Implementation notes:**

- Cubic spline requires uniform nodes — 10 equidistant points were
  manually extracted from the characteristic curve (see report Fig. 1)
- Degree-5 polynomial required input scaling (`u_scaled = u / u_max`)
  to avoid ill-conditioned Vandermonde matrices
- Without clamping M(u) ∈ [0.1, 1.0] H and u₁ ∈ [20, 300] V,
  Lagrange interpolation outside the training range causes simulation
  blow-up (Runge's phenomenon) — a concrete demonstration of
  out-of-distribution failure in data-driven models

All four methods converge to the same solution under the 240·sin(t)
forcing signal when bounds are enforced.

---

## Part 3 — Numerical Integration: Energy Aggregation

Total energy dissipated on resistors R₁ and R₂ over 30 seconds:
W = ∫₀³⁰ [R₁·i₁²(t) + R₂·i₂²(t)] dt

Two composite integration methods, two step sizes (h₁=0.001s, h₂=0.4s),
both linear and nonlinear M models:

**Composite midpoint rule:**
W ≈ h · Σ p((tₖ + tₖ₊₁)/2)
Midpoint values are taken directly from the Improved Euler half-step,
requiring no extra function evaluations.

**Composite Simpson's rule:**
W ≈ (h/3) · [p(t₀) + 4p(t₁) + 2p(t₂) + ... + 4p(t_{m-1}) + p(t_m)]

**Selected results (linear M):**

| Forcing | Midpoint h=0.001 | Simpson h=0.001 | Midpoint h=0.4 | Simpson h=0.4 |
|---------|-----------------|----------------|---------------|--------------|
| 1 V DC | 0.1169 J / 0.0039 W | 0.1169 J / 0.0039 W | 0.1225 J / 0.0041 W | 0.1185 J / 0.0040 W |
| Square wave | 1444.75 J / 48.16 W | 1444.75 J / 48.16 W | 2501.99 J / 83.40 W | 2405.99 J / 80.20 W |
| 240·sin(t) | 184164 J / 6138.8 W | 184164 J / 6138.8 W | 240739 J / 8024.6 W | 231151 J / 7705 W |

At h=0.001 both methods give identical results.
At h=0.4, Simpson's rule consistently outperforms midpoint —
same principle as higher-order aggregation reducing error
in downsampled time-series.

---

## Part 4 — Root-Finding: Frequency Calibration

**Problem:** Find forcing frequency f such that average power P(f) = 406 W.
F(f) = P(f) − 406 = 0

Each evaluation of F(f) runs a full 30-second ODE simulation.
This is the structure of any **expensive black-box optimization**
(e.g. hyperparameter tuning, model calibration).

Three root-finding algorithms benchmarked:

| Method | f found [Hz] | Iterations | F(f) evaluations |
|--------|-------------|-----------|-----------------|
| Bisection | 0.2709 | 10 | 12 |
| Secant | 0.2710 | 8 | 9 |
| Quasi-Newton | 0.2710 | 7 | 14 |

**Quasi-Newton** uses finite-difference derivative approximation:
dF/df ≈ [F(f + Δf) − F(f)] / Δf,   Δf = 1e-4
Step Δf was validated so that halving it changes the estimate by < 1%.
Fewest iterations, but two F(f) calls per step — tradeoff between
iteration count and evaluation cost identical to gradient approximation
in numerical optimization.

**Physical result:** f = 0.271 Hz matches the circuit's natural resonance
frequency, consistent with C=0.5F, L₁=3H, L₂=5H.

---

## How to Run

Requires MATLAB (tested on R2021b+). No additional toolboxes needed.

```matlab
% Run each part independently
czesc1   % ODE solver comparison — generates 8 figures
czesc2   % Nonlinear inductance simulation — generates 8 figures
czesc3   % Energy integration + prints results table to console
czesc4   % Frequency root-finding + prints comparison table to console
```

Each script is self-contained (`clear; close all` at the top).

---

## Skills Demonstrated

- ODE system derivation and state-space formulation
- Explicit numerical integration — stability analysis, step size selection
- Polynomial interpolation (Lagrange, B-spline) and least-squares
  approximation — handling sparse, non-uniform measurement data
- Numerical integration of signals — midpoint and Simpson's rules
- Iterative root-finding with expensive objective functions —
  bisection, secant, quasi-Newton with finite-difference gradients
- Modular MATLAB code structure — reusable `pochodne()`,
  `interpolacja_splajnami()`, `calka_parabol_c3()` functions

---

## Report

Full experiment report with all figures, result tables, and analysis
available in [`Raport.pdf`](Raport.pdf).

---

## License

MIT
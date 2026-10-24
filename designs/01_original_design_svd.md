# Multivariable Optimization of CHO Fed-Batch Parameters — 7 L / 15 L Bioreactors

*Original design (v1) · 2026-10-23 · Ji Chen · Draft proposal and design*

> Live, editable version with figures: [Claude Doc](https://claude.ai/artifact/MEE9WVGSyKkfqxD31KCbKT). The five figures are not in this Markdown copy; their places are marked *[Figure]*.

## Summary

We will set the culture parameters for the new 7 L and 15 L bioreactors with one designed multivariable study, not by one-factor transfer from the shake-flask platform. We treat the process as a map from a scaled input vector **x** (pH, DO, agitation, gas flow, temperature shift, feed rate) to a scaled response vector **y** (ΔVCD, viability, Δlactate, titer, and more). The singular value decomposition (SVD) of the fitted sensitivity matrix shows the few combined input directions that control most of the response ("stiff" modes) and the directions the process does not care about ("sloppy" modes).

Expected outputs:

- A ranked set of 2–3 control modes, each with its input recipe and its response signature.
- A quadratic response model with canonical analysis, which shows if the operating point is a maximum, a saddle, or a ridge.
- A proposed setpoint and NOR/PAR for each vessel, confirmed in one run.
- A reusable JSON data model and Python pipeline for later programs.

## Problem framing

The shake-flask platform is stable, but it cannot tell us the bioreactor setpoints. A flask has surface aeration, no pH or DO control, and a different CO₂ balance. The stirred tank adds new control variables, and these variables interact.

One-factor-at-a-time (OFAT) transfer has three problems:

- It needs many sequential 14-day runs, and we have few vessels.
- It cannot see interactions. Example: a pH setpoint change shifts lactate, base addition, and osmolality, and the effect depends on the DO and gas-flow settings.
- It gives a single point, not a region. It does not tell us how wide the NOR and PAR can be.

A multivariable map answers a different question: *which combinations of inputs move which combinations of outputs, and by how much?* The biology from the flask platform (medium, feed, seeding density, temperature shift) stays as the baseline. Only the physical and control parameters go into the study.

## Mathematical framework

The study models the bioreactor as a map from settings to responses. Near the base operating point, we use a local linear model:

```math
\tilde{y} \approx J\,\tilde{x}, \qquad J_{ij} = \frac{\partial \tilde{y}_i}{\partial \tilde{x}_j}
```

**Step 1 — Scale both spaces.** Inputs and outputs have different units, so we scale them before any decomposition. Without scaling, the singular values only show the choice of units.

```math
\tilde{x}_j = \frac{x_j - c_j}{h_j} \in [-1, +1], \qquad \tilde{y}_i = \frac{y_i - \bar{y}_i}{s_i}
```

cⱼ is the center point and hⱼ the half-range of input j. sᵢ is the SD of output i or, preferred, half of its acceptance or target range. Then one unit of output means "one tolerance", and all outputs have equal weight.

**Step 2 — Decompose J.**

```math
J = U\,\Sigma\,V^{\top} = \sum_{k} \sigma_k\, u_k\, v_k^{\top}, \qquad \sigma_1 \ge \sigma_2 \ge \dots \ge 0
```

*[Figure: SVD of the sensitivity matrix · J = UΣVᵀ]*

Each mode k is one rank-1 piece of J: an input recipe, a gain, and a response signature.

| Object | Size | Meaning in this study |
| --- | --- | --- |
| J | m × p | Sensitivity of each scaled output to each scaled input |
| vₖ (columns of V) | p | Input recipe: a combined setpoint change, e.g. +0.7 pH, −0.5 DO, +0.3 rpm |
| σₖ (diagonal of Σ) | scalar | Gain: how far the responses move per unit of that recipe |
| uₖ (columns of U) | m | Response signature: the pattern of output changes that the recipe causes |

**Why not eigenvalues directly.** J maps a p-dimensional space onto a different m-dimensional space, so Jv = λv is not defined. The eigenvalue problem appears one step removed: JᵀJ = VΣ²Vᵀ. Its eigenvectors are the input recipes, and its eigenvalues are σₖ².

*[Figure: Geometric reading · unit circle to response ellipse]*

The long axis of the ellipse is the direction to control tightly; the short axis is where a wide range costs little.

**Effective dimension.** Mode k carries the share σₖ² / Σσ² of the scaled response variance. The leading modes that together carry ≥ 90% (proposed threshold) form the stiff set. The other modes are sloppy.

## Nonlinear extension: quadratic model and canonical analysis

The linear J is valid only near one operating point, and CHO responses are curved (for example, lactate vs. pH, or viability vs. pCO₂). For each key scaled output, or for a combined desirability D, we fit a second-order model:

```math
\tilde{y} = b_0 + b^{\top}\tilde{x} + \tilde{x}^{\top} B\,\tilde{x}
```

B is a symmetric p × p matrix: the diagonal holds the quadratic terms and the off-diagonal holds half of each interaction term. B maps the input space onto itself, so here the eigenvalue problem is exactly correct:

```math
B = W\,\Lambda\,W^{\top}, \qquad \tilde{x}_s = -\tfrac{1}{2}\,B^{-1} b
```

xₛ is the stationary point. The canonical axes wₖ and eigenvalues λₖ tell the shape of the surface around it.

*[Figure: Canonical analysis · three response-surface shapes]*

The ridge case is common in fed-batch data and it is useful: it lets us choose a setpoint far from a failure edge (for example, high osmolality or high pCO₂) at almost no loss of response.

| Eigenvalues of B | Shape | Action |
| --- | --- | --- |
| All λ < 0 | Maximum | xₛ is the optimum, if it is inside the design region |
| All λ > 0 | Minimum | Move away; the optimum is on the region boundary |
| Mixed signs | Saddle | Ridge analysis along the rising direction |
| One or more λ ≈ 0 | Ridge | Pick the most robust point on the ridge |

**Link to the SVD.** The local sensitivity at any point is the gradient, J(x̃) = b + 2Bx̃ for each output. After we choose the setpoint, we repeat the SVD of J at that point. The NOR/PAR is set from that local SVD, not from the one at the design center.

## Experimental design

We use a definitive screening design (DSD) with vessel scale as a 2-level categorical factor. A DSD tests each continuous factor at three levels in about 2k+1 runs. Main effects are not aliased with two-factor interactions, and curvature is detectable. For 7 factors this is about 17 runs; adding 2 "fake" factors (about 21 runs) improves the error estimate.

**Factors (proposed; center = current platform value where one exists)**

| Factor | Unit | Low / center / high | Basis |
| --- | --- | --- | --- |
| pH setpoint | – | 6.85 / 7.00 / 7.15 | Common CHO range; confirm against platform |
| DO setpoint | % air sat. | 30 / 45 / 60 | Common CHO range |
| Agitation, as P/V | W/m³ | From vessel characterization | Low = minimum for peak O₂ demand; high = shear or foam limit |
| Sparge rate | vvm | From vessel characterization | Low = minimum for CO₂ stripping; high = foam limit |
| Temperature-shift day | day | Platform − 1 / platform / platform + 1 | Or a VCD trigger, if the platform uses one |
| Feed volume | % of platform | 80 / 100 / 120 | Same feed schedule, scaled bolus volume |
| Vessel | categorical | 7 L / 15 L | Tests if the setpoints transfer between the two systems |

**Fixed factors:** medium and feed formulation, seed train, seeding density, glucose target, antifoam strategy, harvest day, and base type. These come from the flask platform and stay constant.

**Blocking.** Runs go in blocks set by the number of vessels available at the same time. Each block uses one seed train and includes one center-point run. The center points measure block-to-block drift and pure error.

**Stages**

1. Screening DSD: finds the active factors and first curvature.
2. Augmentation: adds runs around the best region only for the 2–4 active factors (or uses Bayesian optimization when vessels are few).
3. Confirmation: one run per vessel at the chosen setpoint.

## Data model

Each run produces one input record and one output record, joined by `run_id`. Records hold the **actual** achieved values, not only the setpoints, because the model must use what the culture saw.

**Input record**

```json
{
  "run_id": "R07",
  "block": 2,
  "vessel": { "id": "BR15-01", "scale_L": 15 },
  "seed_train": "ST-2026-11",
  "factors_setpoint": {
    "pH_sp": 7.00, "DO_sp": 45, "P_V_Wm3": null,
    "sparge_vvm": null, "T_shift_day": 5, "feed_pct": 100
  },
  "factors_actual": {
    "pH_mean": null, "DO_mean": null, "P_V_Wm3": null,
    "sparge_vvm_mean": null, "T_shift_day": null, "feed_pct": null
  },
  "fixed": { "medium": "", "feed": "", "seed_VCD_e6": null, "harvest_day": 14 }
}
```

**Output record** (batch-level summaries for the model; the daily series is kept for later work)

```json
{
  "run_id": "R07",
  "responses": {
    "dVCD_growth_e6": null, "IVCD_e9_cell_day_L": null,
    "viability_harvest_pct": null, "dLactate_mM": null,
    "titer_gL": null, "qP_pg_cell_day": null,
    "osmo_max_mOsm": null, "pCO2_max_mmHg": null, "NH4_max_mM": null,
    "HMW_pct": null, "main_peak_pct": null, "afucosylation_pct": null
  },
  "timeseries": [ { "day": 0, "VCD_e6": null, "viab_pct": null, "glc_gL": null, "lac_mM": null } ]
}
```

**Response definitions.** Each response has one fixed definition, a goal, and a scale sᵢ for the scaling step.

| Response | Definition | Goal | Scale sᵢ |
| --- | --- | --- | --- |
| ΔVCD (growth) | VCD at temperature shift − VCD at day 0 | Maximize | Platform SD or target window |
| Viability at harvest | % viable on harvest day | Maximize, ≥ limit | Half of acceptance range |
| ΔLactate | Lactate at harvest − peak lactate (negative = consumption) | Minimize | Platform SD |
| Titer | g/L at harvest | Maximize | Platform SD |
| Osmolality, max | Highest value in the run | Constraint (≤ limit) | Distance to limit |
| pCO₂, max | Highest offline value | Constraint (≤ limit) | Distance to limit |
| HMW, main peak, afucosylation | At harvest | Target window | Half of target window |

The goals and limits are open items for the team to set before the first run.

## Analysis pipeline

The same pipeline runs after each block of runs, so the model improves as data arrive.

*[Figure: Analysis pipeline · 7 steps, 1 loop — 1 JSON records → 2 Validate, derive → 3 Scale → 4 Fit models → 5 Decompose → 6 Decide → 7 Confirm → Lock parameters; a failed confirmation → Augment design → new runs → step 1]*

A failed confirmation does not restart the study; it adds runs where the model was wrong.

- **Fit.** One model per output, with linear and quadratic terms chosen by forward selection under effect heredity. PLS2 on all outputs is the cross-check: its weight vectors should point near the leading vₖ.
- **Uncertainty.** Bootstrap over whole runs, stratified by block. This gives confidence intervals for each σₖ, for the loadings of each vₖ, and for xₛ.
- **Validation.** Leave-one-block-out prediction (Q²) and a lack-of-fit test from the center points.

## Decision rules

The model drives a decision only when it passes each gate below. All thresholds are proposed and open for team review.

1. **Model is fit for use.** Leave-one-block-out Q² ≥ 0.5 for each decision response, and no significant lack of fit.
2. **A mode is real.** The bootstrap 95% CI of σₖ is above the noise floor from the replicate center points. Otherwise the mode is noise.
3. **Stiff set.** The smallest set of real modes that carries ≥ 90% of Σσ².
4. **Setpoint.** Maximize the overall desirability D (weighted geometric mean of the per-response dᵢ), with hard limits on osmolality, pCO₂, and viability, inside the design region. On a ridge, take the point farthest from the limits.
5. **PAR and NOR.** From the local J at the chosen setpoint (rule below).
6. **Vessel transfer.** If the vessel factor or its interactions are significant, set vessel-specific setpoints with the same method. Otherwise use one recipe for both.
7. **Confirmation passes** when each decision response falls inside its 95% prediction interval in both vessels.

**PAR rule.** Outputs are scaled in tolerance units, so a coded move δ changes the responses by Jδ. The PAR is the set of moves that keeps the combined change within one tolerance:

```math
\lVert J\,\delta \rVert_{2} \le 1 \;\Longleftrightarrow\; \delta^{\top} V\,\Sigma^{2}\,V^{\top}\,\delta \le 1
```

The half-width along vₖ is 1/σₖ. The NOR uses κ/σₖ, with κ = 0.5 proposed. The 2-norm is conservative: the exact per-output rule, max |(Jδ)ᵢ| ≤ 1, gives a larger polytope that contains this ellipse.

*[Figure: PAR in input space · ellipse, design cap, factor box]*

Operators set factors one at a time, so the reported ranges are the largest upright box inside the ellipse and the design region. When a stiff mode loads on two factors together (for example pH and DO), the box gets thin. Then we either tighten control on that mode or state the coupled limit as an explicit rule.

## Risks and limitations

The main risk is a model that looks clean but is built on too few, correlated runs.

| Risk | Effect | Mitigation |
| --- | --- | --- |
| Few runs for many factors | Overfit; unstable vₖ | DSD, effect sparsity, bootstrap CIs, leave-one-block-out validation |
| Linear J is local | Wrong modes far from the center | Quadratic model; repeat the SVD at the chosen setpoint |
| Control loops couple inputs (DO cascade changes gas flow; pH uses CO₂ and base) | Realized inputs are correlated | Model on actual values; record gas and base totals as covariates |
| Block effects (seed train, media lot) | Drift read as a factor effect | Blocks with center points; block term in the model |
| Output scaling sᵢ is a judgment | Ranking of modes depends on it | Agree sᵢ before the first run; test the sensitivity of the conclusions |
| Batch summaries lose time structure | Miss phase-specific effects (e.g. lactate shift timing) | Keep daily series; phase-wise models later |
| 7 L and 15 L are close in scale | Results do not extrapolate to large scale | Re-use the framework at scale-up with physical features (P/V, kLa, v_s) |
| Product-quality assays measured once, at harvest | Assay noise can mask small effects | Use assay precision in sᵢ; replicate at center points |

## Deliverables, timeline, and next step

The study length depends mostly on how many vessels can run at the same time; each block takes about 3 weeks (14-day run plus turnaround).

| Phase | Work | Estimate |
| --- | --- | --- |
| 0 | Vessel characterization (kLa, P/V, CO₂ stripping) to fix the agitation and sparge ranges | 1–2 weeks |
| 1 | Design and data-model sign-off: factors, ranges, sᵢ, goals, limits | 1 week, in parallel with phase 0 |
| 2 | Screening DSD, in blocks | About 21 runs ÷ vessels per block × 3 weeks |
| 3 | Analysis: fit, SVD, canonical analysis, setpoint, PAR/NOR | 1 week after the last block |
| 4 | Augmentation (only if needed) and confirmation runs | 3–6 weeks |

**Deliverables:** the design table, the JSON schema, the Python pipeline, a mode report (vₖ, σₖ, uₖ with CIs), the setpoint and PAR/NOR per vessel, and the confirmation report.

**Open questions**

- [ ] How many 7 L and 15 L vessels can run at the same time?
- [ ] Current platform values for the center points (pH, DO, temperature-shift day, feed volume)
- [ ] Goals, limits, and sᵢ for each response, including product-quality targets
- [ ] Is product quality part of the optimization objective, or only a constraint?

**Next step: Python implementation.** Modules: schema validation of the JSON records; DSD generation with blocks; per-output model fitting and PLS cross-check; SVD of J and eigen-decomposition of B with bootstrap CIs; desirability optimization under constraints; PAR/NOR box calculation; and a report page. A synthetic data set with known modes will test the pipeline before real runs arrive.

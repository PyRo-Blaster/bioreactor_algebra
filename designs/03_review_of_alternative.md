# Review of the alternative design — mathematics

*2026-10-24 · Review notes on [02_alternative_design_bayesian.md](02_alternative_design_bayesian.md), with defects it shares with [01_original_design_svd.md](01_original_design_svd.md). Numbers checked with [`scripts/dsd_check.py`](../scripts/dsd_check.py).*

## Short answer

The direction of the alternative is correct: shrinkage replaces stepwise selection, decisions use the distance to the limits, and the NOR comes from measured control error. But its mathematics has four real defects, and it keeps two defects from the original design.

## What is sound

- **Shrinkage in place of forward selection.** With this ratio of runs to terms, stepwise selection is unstable. A horseshoe prior gives the same sparsity assumption without a hard cut-off.
- **Draw-by-draw propagation** of uncertainty through the SVD, the stationary point, and D.
- **Chance constraints (Eq. 25) and the NOR (Eq. 27)** use the distance to the specification limits. The original PAR (‖Jδ‖ ≤ 1) assumes the setpoint sits at the center of each tolerance, which is usually false.
- **Measured control error in place of κ = 0.5.**

## Defects, most important first

**1. The 13-run DSD is saturated, so the prior sets the curvature.** The model matrix is 13 × 28 with rank 13, so a 15-dimensional null space is set only by the horseshoe prior. The risk of trusting the prior moves from the design to the analysis.

- Shrinkage pulls B and its eigenvalues toward zero, so the prior creates "ridges" by itself.
- When B⁽ˢ⁾ is near-singular, x̃ₛ = −½B⁻¹b becomes very large, and its posterior has heavy tails. Report the constrained optimum inside the design region, or a ridge analysis.
- Eq. 19 has independent λⱼ, so it does not give effect heredity. Couple λ for each interaction to λ for its parent main effects.
- The fit has no error degrees of freedom. Replicated center points are extra runs, and neither document counts them. With them, the run saving drops to about 10–15%.

**2. The bridge sits where it measures least.** At the corner (+1,…,+1), 45% of φ is outside the space the 13 runs can estimate. The vessel contrast there is confounded with the main vessel's prediction error. Three runs give only 3 of the 7 vessel terms, so a 4-dimensional subspace is invisible. At equal vvm, v_s in the 15 L vessel is about 1.29× that in the 7 L vessel, so a vessel × sparge interaction is built in.
*Fix:* code the gas factor as v_s (or match on kLa or CO₂ stripping), and repeat existing DSD runs in the second vessel. Each pair is then a model-free paired difference with variance 2σ².

**3. The pass-probability gate gives the same result every time.** If the criterion is "inside its own 95% predictive interval," then P(all pass) ≈ 0.95¹² = 0.54 for 6 responses × 2 vessels. That is always in the 0.50–0.80 band. The original confirmation rule has the same multiplicity problem: a perfect model fails about 46% of the time.
*Fix:* gate on P(y_new ∈ spec). Check the model jointly, with a Mahalanobis distance against χ² or with 99.6% per-interval levels (Bonferroni).

**4. CVaR combines epistemic and aleatory uncertainty.** D(x; θ) varies only through θ, so CVaR favors regions where there is more data, not regions where the process is flat. Epistemic uncertainty should decide where to *sample*, and aleatory uncertainty (control error, batch noise) should decide where to *operate*. If any dᵢ can be 0, CVaR₀.₁₀ is flat at 0 and the argmax has no unique solution.

**5. The PAR in Eq. 26 is not a coherent bound.** Sorted singular values swap modes between draws, the marginal bounds σₖ,up are not simultaneous, and the uncertainty in V̂ is ignored.
*Cleaner:* M = E[JᵀJ | Y] = J̄ᵀJ̄ + Σᵢ Cov(Jᵢ·), and PAR: δᵀMδ ≤ 1. For the modes, report the stiff subspace (projector V_K V_Kᵀ and principal angles), not individual vectors. Apply the distance-to-limit idea to the PAR too.

**6. Shared with the original: gas and base totals are mediators, not covariates.** The path pH setpoint → lactate → base → osmolality means that conditioning on base removes part of the pH effect. Treat these totals as responses or in a path model. Related: pH probe drift of about ±0.05 is 0.33 coded units, which attenuates the pH coefficients. Use offline-corrected pH.

**Minor points.** Section 5 gives ARD length scales per factor, but Eq. 21 gives them per stiff coordinate. z is derived from the same data the GP is fitted to. "Largest R" in Eq. 27 needs a shape family (for example a box scaled by t, with t maximized). A block random effect from 3–5 blocks is mostly set by its prior.

## Verdict

| Item | Decision |
|---|---|
| Horseshoe shrinkage, draw-by-draw propagation | Keep. Add heredity coupling and prior-sensitivity checks (τ₀, c) |
| Chance constraints; NOR from measured control error | Keep. Apply distance to limits to the PAR |
| 13-run DSD with no error df | Fix. Count the center runs, or go back to fake factors |
| Bridge at diagonal corners | Replace with replicated DSD runs; code gas as v_s |
| Pass-probability gate; confirmation rule | Fix multiplicity (both designs) |
| CVaR over θ only | Change: aleatory risk for the setpoint, epistemic uncertainty for sampling |
| Gas and base totals as covariates | Remove from the regressors (both designs) |

Recommended path: put the alternative's analysis and decision rules on top of a screening design with enough runs, with the corrections above.

## Interdisciplinary note

CVaR comes from finance, where Rockafellar and Uryasev used it to hedge market risk, which you cannot learn away. Applying it to epistemic uncertainty is like insuring against your own ignorance when information costs only a few runs. Reinforcement learning uses opposite signs: epistemic uncertainty gets an optimism bonus (UCB exploration), and aleatory uncertainty gets a pessimism penalty. A sequential study could change sign over time: optimistic in Round 2 while it learns, pessimistic in Round 3 when it commits. Belevitch studied conference matrices, the basis of every DSD, around 1950 for lossless conference-telephone networks, where every line hears every other line equally. That is the property that keeps the main effects orthogonal. The bridge problem resembles ICH E5 bridging: bridge with a matched design at points you already know, not at the extremes.

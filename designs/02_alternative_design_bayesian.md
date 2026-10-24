**Multivariable Optimization of CHO Fed-Batch Parameters**

Alternative Design — Fewer Runs, Robustness-First • 7 L / 15 L Bioreactors

*Draft alternative proposal for side-by-side review with the original design • October 2026*

How to read this document. This is a complete alternative, not a set of edits. It keeps the original study's question, factor set, scaling, SVD / canonical-analysis framework, actual-value data model and JSON pipeline. It changes three things only: the size and structure of the screening design, the analysis (shrinkage and hierarchical Bayesian fitting instead of stepwise selection), and the decision rules (robust, posterior-based setpoint, PAR and NOR). Section 2 sets out the mathematical framework in full — the scaling, the matrices J and B, the SVD and canonical analysis from the original's theory section, written out as matrices and formulas, plus the posterior extensions this design adds. Section 9 gives the head-to-head run and timeline comparison.

# Summary — the same questions in about 21 runs and 3 culture rounds

The original design answers the right questions but pays for them twice: once with a 21-run screening design that crosses vessel as a full factor, and again with open-ended augmentation. This alternative keeps a classical, model-lean screening skeleton — so the design itself does not depend on the priors being right — and moves the Bayesian machinery into analysis and decisions, where it reduces run count without reducing robustness.

- **Round 1 — Screening (16 runs).** 13 runs in the main vessel (a standard definitive screening design, DSD, for 6 continuous factors) plus 3 bridge runs in the second vessel — 16 runs in Round 1, instead of about 21.

- **Round 2 — Adaptive correction (3–4 runs).** 3–4 runs placed only where the decision is sensitive or two models disagree, gated by a pre-computed probability of passing confirmation. This round can shrink to 2 runs, or merge into confirmation, if that probability is already high.

- **Round 3 — Confirmation (2 runs).** 1 run per vessel at the chosen setpoint, judged against pre-registered 95% posterior predictive intervals.

Total: about 21 culture runs in 3 culture rounds, compared with 25–30 runs once open-ended augmentation is included in the original plan. The saving is deliberately modest. Six continuous factors have an information floor: pushing Round 1 down to 11–12 runs with a purely Bayesian optimal design would trade robustness for count, which is the opposite of the brief for this alternative.

| *Design principle in one line: keep the design classical and robust; make the analysis and the decisions Bayesian. Priors inform the fit and the choice of setpoint — they do not determine where the runs go in Round 1.* |
|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|

# 1 Design principles — where Bayesian methods help, and where they do not

1.  **Do not gamble the design on the prior.** A purely Bayesian D-optimal design is efficient only if the prior is right. Platform experience in shake flasks may not transfer cleanly to stirred tanks, so Round 1 does not bet the run placement on it.

2.  **Shrink, do not select-and-discard.** With roughly 21 runs and 35 candidate terms in a full quadratic model (7 linear terms counting vessel handling, 6–7 pure quadratics and 15–21 interactions, depending on coding), stepwise forward selection is unstable: a different batch can select a different term set. Shrinkage estimation keeps all plausible terms and pulls unsupported ones toward zero.

3.  **Decide under uncertainty, not on a point estimate.** The setpoint, PAR and NOR are chosen from the full posterior, not from a point estimate. A flat ridge that looks optimal on a point estimate can be fragile once parameter uncertainty is included.

4.  **Pre-register everything that can be gamed later.** Decision rules, stopping rules and confirmation criteria are written and agreed before Round 1 starts. Thresholds are not moved after data arrive.

What does not change: inputs and outputs are scaled exactly as in the original (input coded to \[−1, +1\]; output scaled so one unit equals one tolerance, sᵢ = half the acceptance or target range where available); models are fitted on actual achieved values, not setpoints; gas and base totals are recorded as covariates; and the local Jacobian J at the chosen setpoint — not at the design centre — drives PAR and NOR.

# 2 Mathematical framework — the shared map, its matrices, and the posterior extension

This section reproduces the original proposal's theoretical explanation in full — the same scaling, the same sensitivity matrix J, the same singular value decomposition (SVD) and the same canonical analysis of the quadratic surface — because this alternative changes none of that mathematics. What it adds is stated explicitly at the end (Section 2.4): in the original, J and B are point estimates; here they are posterior distributions, and every quantity derived from them — the gains σₖ, the stationary point, the desirability, PAR and NOR — inherits that distribution. Equation numbers are used for cross-reference in Sections 5–8.

## Scaling, and the local linear map

The study models the bioreactor as a map from settings to responses. Step 1 is to scale both spaces. Inputs and outputs have different units, so they are scaled before any decomposition; without scaling, the singular values would only show the choice of units. Input j is coded by its centre cⱼ and half-range hⱼ, and output i by sᵢ — the standard deviation of output i or, preferred, half of its acceptance or target range. Then one unit of scaled output means “one tolerance”, and all outputs carry equal weight:


```math
\tilde{x}_j = \frac{x_j - c_j}{h_j} \in [-1,\,+1], \qquad \tilde{y}_i = \frac{y_i - \bar{y}_i}{s_i} \tag{1}
```


Step 2 is the local linear model. Near the base operating point the map is approximated by its Jacobian J — an m × p matrix (m scaled responses, p = 6 continuous factors in this study), whose entry Jᵢⱼ is the sensitivity of scaled output i to scaled input j. Written out in full, J is the matrix of partial derivatives:


```math
\tilde{y} \approx J\,\tilde{x}, \qquad J_{ij} = \frac{\partial \tilde{y}_i}{\partial \tilde{x}_j} \tag{2}
```


*(Equation 3: J written out as a matrix — see the source .docx.)*

| **Object**         | **Size** | **Meaning in this study**                                                       |
|--------------------|----------|---------------------------------------------------------------------------------|
| J                  | m × p    | Sensitivity of each scaled output to each scaled input                          |
| vₖ (columns of V)  | p        | Input recipe: a combined setpoint change, e.g. +0.7 pH, −0.5 DO, +0.3 agitation |
| σₖ (diagonal of Σ) | scalar   | Gain: how far the responses move per unit of that recipe                        |
| uₖ (columns of U)  | m        | Response signature: the pattern of output changes that the recipe causes        |

## Decomposing J — the SVD, stiff modes and sloppy modes

J is decomposed by the SVD, with Σ the diagonal matrix of gains ordered from largest to smallest:


```math
J = U\,\Sigma\,V^{\top} \tag{4}
```


*(Equation 5: \Sigma written out as a matrix — see the source .docx.)*

Equivalently, J is a sum of rank-1 pieces — one per mode. Each mode k is an input recipe vₖ, a gain σₖ, and a response signature uₖ: a unit move along the recipe moves the responses by exactly the gain, in the signature's direction, and a general move δ decomposes over the recipes:


```math
J = \sum_{k=1}^{r} \sigma_k\, u_k\, v_k^{\top}, \qquad \sigma_1 \geq \sigma_2 \geq \cdots \geq \sigma_r &gt; 0 \tag{6}
```



```math
J v_k = \sigma_k u_k, \qquad \delta = \sum_k a_k v_k \;\Longrightarrow\; J\delta = \sum_k a_k\, \sigma_k\, u_k \tag{7}
```


Why not eigenvalues of J directly? J maps a p-dimensional space onto a different m-dimensional space, so Jv = λv is not defined. The eigenvalue problem appears one step removed, on the two Gram matrices, whose eigenvectors are exactly the recipes and signatures and whose eigenvalues are the squared gains:


```math
J^{\top} J = V\,\Sigma^{2}\,V^{\top}, \qquad J J^{\top} = U\,\Sigma^{2}\,U^{\top} \tag{8}
```


Geometric reading: J maps the unit ball of coded input moves onto an ellipsoid in response space whose semi-axes have lengths σₖ and point along uₖ. The long axis is the direction to control tightly; the short axis is where a wide range costs little. Mode k carries the share ρₖ of the scaled response variance, and the stiff set is the smallest set of leading modes that together carries at least 90% (proposed threshold, as in the original); the remaining modes are sloppy:


```math
\rho_k = \frac{\sigma_k^{2}}{\sum_j \sigma_j^{2}}, \qquad \mathrm{stiff\ set:\ smallest\ } K \mathrm{\ with\ } \sum_{k=1}^{K} \rho_k \geq 0.90 \tag{9}
```


## Nonlinear extension — the quadratic model and canonical analysis

The linear J is valid only near one operating point, and CHO responses are curved (for example, lactate vs. pH, or viability vs. pCO₂). For each key scaled output, or for a combined desirability D, a second-order model is fitted over the coded factors. B is a symmetric p × p matrix: the diagonal holds the pure quadratic coefficients and the off-diagonal holds half of each interaction coefficient, because each interaction appears twice in the quadratic form:


```math
\tilde{y} = b_0 + b^{\top}\tilde{x} + \tilde{x}^{\top} B\,\tilde{x} \tag{10}
```


*(Equation 11: B written out as a matrix — see the source .docx.)*

The candidate term count follows directly from the matrix: one intercept, p linear terms, p pure quadratics on the diagonal of B, and one term per off-diagonal pair. For p = 6 that is 28 terms per response before any vessel terms are added — the count that drives the analysis choice in Section 5:


```math
\mathrm{number\ of\ terms} = 1 + p + p + \frac{p(p-1)}{2} = 28 \quad (p = 6) \tag{12}
```


Unlike J, B maps the input space onto itself, so here the eigenvalue problem is exactly correct. Differentiating (10) gives the gradient and the Hessian, and setting the gradient to zero gives the stationary point x̃ₛ. The local sensitivity at any point is that gradient — for a vector of responses, the rows of J(x̃) — so after a setpoint is chosen, the SVD of Section 2.2 is repeated at that point, not at the design centre. That local J, not the centre J, sets PAR and NOR (Section 6):


```math
\nabla \tilde{y}(\tilde{x}) = b + 2B\tilde{x}, \qquad \nabla^{2}\tilde{y} = 2B \tag{13}
```



```math
\tilde{x}_s = -\frac{1}{2}\, B^{-1} b, \qquad B = W\,\Lambda\,W^{\top} \tag{14}
```


In canonical coordinates z = Wᵀ(x̃ − x̃ₛ) the surface separates into independent one-dimensional parabolas, one per eigenvalue λₖ of B. The signs of the λₖ classify the surface, exactly as in the original:


```math
\tilde{y} = \tilde{y}_s + \sum_{k=1}^{p} \lambda_k z_k^{2}, \qquad z = W^{\top}(\tilde{x} - \tilde{x}_s) \tag{15}
```


| **Eigenvalues of B** | **Shape** | **Action**                                           |
|----------------------|-----------|------------------------------------------------------|
| All λ \< 0           | Maximum   | x̃ₛ is the optimum, if it is inside the design region |
| All λ \> 0           | Minimum   | Move away; the optimum is on the region boundary     |
| Mixed signs          | Saddle    | Ridge analysis along the rising direction            |
| One or more λ ≈ 0    | Ridge     | Pick the most robust point on the ridge              |

The ridge case is common in fed-batch data and it is useful: it lets the team choose a setpoint far from a failure edge (for example, high osmolality or high pCO₂) at almost no loss of response. Section 6 makes that choice quantitative instead of visual.

## The Bayesian extension — the same matrices, as distributions

Everything above is kept. The change is that (b, B) are no longer single fitted values. Write the quadratic model for response i and run r in regression form, with φ(x̃) the vector of the 28 candidate terms in (12), a block (seed-batch) random effect γ, and residual error ε. Stacking the runs gives the design matrix X whose rows are φ(x̃ᵣ)ᵀ:


```math
y_{ir} = \varphi(\tilde{x}_r)^{\top} \beta_i + \gamma_{b(r),\,i} + \varepsilon_{ir}, \qquad \gamma_{bi} \sim \mathcal{N}(0,\, \tau_b^{2}), \quad \varepsilon_{ir} \sim \mathcal{N}(0,\, \sigma_i^{2}) \tag{16}
```



```math
\varphi(\tilde{x}) = \left[\, 1,\ \tilde{x}_1, \ldots, \tilde{x}_p,\ \tilde{x}_1^{2}, \ldots, \tilde{x}_p^{2},\ \tilde{x}_1\tilde{x}_2, \ldots, \tilde{x}_{p-1}\tilde{x}_p \,\right]^{\top} \tag{17}
```


Bayes' rule then turns the likelihood of the run data Y and the priors registered in Round 0 (Section 3) into a posterior over every coefficient at once. There is no selected term set: the regularized horseshoe prior below keeps supported effects large and shrinks unsupported ones smoothly toward zero, with the shrinkage of an interaction tied to its parent main effects so that effect heredity holds:


```math
p(\beta, \gamma, \sigma, \tau \,|\, Y) \;\propto\; p(Y \,|\, \beta, \gamma, \sigma)\; p(\beta \,|\, \tau, \lambda)\; p(\gamma \,|\, \tau_b)\; p(\tau, \lambda, \sigma) \tag{18}
```



```math
\beta_j \,|\, \lambda_j, \tau \sim \mathcal{N}\!\left(0,\; \tau^{2} \tilde{\lambda}_j^{2}\right), \qquad \tilde{\lambda}_j^{2} = \frac{c^{2}\lambda_j^{2}}{c^{2} + \tau^{2}\lambda_j^{2}}, \qquad \lambda_j,\, \tau \sim \mathrm{half\!-\!Cauchy} \tag{19}
```


Every derived object is computed draw by draw. For posterior draw s, the draw's (b⁽ˢ⁾, B⁽ˢ⁾) gives a draw of the local Jacobian at the candidate setpoint x̃\*, of its SVD, and of the stationary point — so σₖ, vₖ, uₖ and x̃ₛ are reported as posterior medians with credible intervals, not as single numbers with a separate bootstrap:


```math
J^{(s)}(\tilde{x}^{*}) = U^{(s)}\,\Sigma^{(s)}\,V^{(s)\top}, \qquad \tilde{x}_s^{(s)} = -\frac{1}{2}\left(B^{(s)}\right)^{-1} b^{(s)}, \qquad s = 1, \ldots, S \tag{20}
```


A Gaussian-process (GP) cross-check is fitted in the stiff subspace only. With z = V_Kᵀx̃ the coordinates along the K stiff recipes from (9), the GP covariance uses one length scale ℓ_d per stiff coordinate (automatic relevance determination, ARD) — a short ℓ_d flags a coordinate the quadratic surface is smoothing over:


```math
z = V_K^{\top}\tilde{x}, \qquad k(z, z&#39;) = \sigma_f^{2}\, \exp\!\left(-\frac{1}{2}\sum_d \frac{(z_d - z_d&#39;)^{2}}{\ell_d^{2}}\right) \tag{21}
```


Finally, decisions are functionals of the posterior, defined here and applied in Sections 6–8. The overall desirability D is the weighted geometric mean of the per-response desirabilities dᵢ, as in the original — but it is now a random variable through the posterior draw θ⁽ˢ⁾ = (b⁽ˢ⁾, B⁽ˢ⁾, σ⁽ˢ⁾). Predictions for a future run integrate over that posterior (the posterior predictive), which is what confirmation is judged against in Section 8:


```math
D(\tilde{x};\, \theta) = \left(\prod_i d_i^{\,w_i}\right)^{\!1/\sum_i w_i} \tag{22}
```



```math
p(\tilde{y}_{\mathrm{new}} \,|\, Y) = \int p(\tilde{y}_{\mathrm{new}} \,|\, \theta)\; p(\theta \,|\, Y)\; d\theta \tag{23}
```


Reading guide for the rest of the document: Sections 3–4 need only the classical objects (1)–(15) to state the design. Sections 5–8 use the posterior objects (16)–(23). Where a formula in the original and a formula here appear to conflict — for example the NOR constant κ in Section 6 — the posterior version is the operative one, and the difference is flagged at the point of use.

# 3 Round 0 — Prior registration (no culture runs)

Before Round 1, the team holds one structured elicitation session (half a day is usually enough) and records, for each continuous factor, the expected direction and plausible magnitude of its first-order effect on each decision response, plus a plausible range for run-to-run and block-to-block standard deviations. Sources are named explicitly: shake-flask platform data, prior reactor runs, and literature values. Disagreement is recorded as a wider prior, not averaged away.

- Priors on first-order effects and noise, in scaled (tolerance) units, so they can be used directly in the shrinkage fit in Section 5.

- The response goals, hard limits and sᵢ values — unchanged from the original, and agreed before the first run, because the scaling determines the stiff / sloppy ranking.

- The thresholds in Sections 6–8 (chance-constraint levels, CVaR level, confirmation-pass probability and interval). All are proposals for team review; once agreed, they are frozen.

If the team cannot agree on quantitative priors, the design still works: the shrinkage fit in Section 5 falls back to weakly informative priors and behaves like regularized regression. The gain from Round 0 is smaller in that case, but nothing in Rounds 1–3 breaks.

# 4 Round 1 — Screening in 16 runs: a 13-run DSD plus a 3-run vessel bridge

Six continuous factors are screened in the main vessel with a standard DSD: 2 × 6 + 1 = 13 runs. A DSD gives three levels per factor, main effects unaliased with two-factor interactions, and a first read on curvature. The original design's additional runs for 'fake' factors are dropped; pure error is estimated from replicated centre points instead, which also serve as the drift monitors in every block.

| **Factor**            | **Unit**      | **Low / centre / high**              | **Notes**                                        |
|-----------------------|---------------|--------------------------------------|--------------------------------------------------|
| pH setpoint           | —             | 6.85 / 7.00 / 7.15                   | As in the original; confirm against platform     |
| DO setpoint           | % air sat.    | 30 / 45 / 60                         | As in the original                               |
| Agitation, as P/V     | W/m³          | From vessel characterization         | Low: O₂-demand minimum; high: shear / foam limit |
| Sparge rate           | vvm           | From vessel characterization         | Low: CO₂-stripping minimum; high: foam limit     |
| Temperature-shift day | day           | Platform −1 / platform / platform +1 | Or a VCD trigger if the platform uses one        |
| Feed volume           | % of platform | 80 / 100 / 120                       | Same schedule, scaled bolus volume               |

Fixed factors are unchanged from the original: medium and feed formulation, seed train, seeding density, glucose target, antifoam strategy, harvest day and base type.

## Vessel is bridged, not crossed

Vessel does not enter the DSD as a fully crossed factor. Instead, the second vessel runs a 3-run bridge: the centre point plus two diagonal extreme points chosen to be maximally informative about a vessel main effect and vessel-by-factor interaction. The 7 L and 15 L systems are geometrically close, so a large interaction is a priori unlikely; the bridge tests that assumption cheaply instead of paying for it in every run.

- If the bridge shows a small vessel effect (posterior interval for the vessel terms includes zero and is narrow relative to one tolerance unit), one setpoint serves both vessels, with the bridge runs pooled into the fit.

- If the bridge shows a large or uncertain vessel effect, Round 2 is re-planned as vessel-specific follow-up (see Section 7) — the cost is paid only when the data say it is needed.

## Each run yields more than its endpoint

In addition to the batch-level responses in the original data model, each run contributes 2–3 pre-specified time-course features extracted from the daily series — for example exponential-phase growth rate and the day of peak lactate. These cost no extra runs and materially improve the precision of the growth and lactate models. The daily series itself is retained in full, as in the original.

Blocking is unchanged in spirit: runs are blocked by simultaneous vessel availability, each block uses one seed train and includes a centre point. Block (seed batch) is modelled explicitly as a random effect (Section 5), so drift becomes a quantified variance component rather than unexplained noise.

# 5 Analysis — shrinkage, hierarchy and two models in parallel

## Shrinkage fitting replaces stepwise selection

Each decision response is fitted with a Bayesian quadratic model over the candidate term set, using a regularized horseshoe prior (with the Round 0 priors where they exist). Effects supported by the data stay large; unsupported effects shrink smoothly toward zero. Effect heredity is respected by construction of the prior: an interaction is allowed a large coefficient only if its parent main effects are not shrunk to zero. This replaces forward selection, whose term set is unstable at this run-to-parameter ratio, and it yields a full posterior for every coefficient — which Sections 6 and 7 need.

## Hierarchy, actual values and covariates

- Block / seed batch enters as a random effect, so block-to-block variance is estimated and propagated into prediction intervals instead of inflating residual error.

- All models use actual achieved factor values. Total gas flow and total base addition per run are included as covariates, because the DO cascade and pH control loops couple the realized inputs (a risk already identified in the original).

- Centre points in every block continue to provide pure error, drift monitoring and a lack-of-fit check.

## Two surface models, and disagreement as information

Two models are carried in parallel: (a) the quadratic response surface above, with canonical analysis of B (eigenvalues, stationary point xₛ = −½ B⁻¹b, ridge / saddle / maximum classification) exactly as in the original; and (b) a Gaussian-process (GP) model fitted in the low-dimensional stiff subspace identified by the SVD of J, with automatic relevance determination (ARD) length scales per factor. The GP is not a replacement for the quadratic model at this sample size — it is a cross-check. Where the two models' predictions diverge materially, that location is informative, and Round 2 runs are placed there first (Section 7). PLS2 on all outputs, as in the original, remains a second cross-check on the leading SVD modes.

Uncertainty is reported from the posterior directly (credible intervals for σₖ, for the loadings of each vₖ, and for xₛ), with leave-one-block-out predictive performance (the analogue of Q²) as the model-adequacy gate. Bootstrap over whole runs, as in the original, is retained as a non-Bayesian cross-check on the headline intervals.

# 6 Robust setpoint, PAR and NOR — decisions that survive a wrong model

## Setpoint: maximize robust desirability, subject to chance constraints

Let D be the overall desirability (weighted geometric mean of the per-response desirabilities dᵢ), as in the original. The setpoint does not maximize D at the point estimate. It maximizes the conditional value at risk (CVaR) of D under the posterior — the mean of D over its worst 10% of posterior draws (level proposed; team to confirm) — subject to chance constraints on the hard limits, each evaluated under the same posterior:


```math
\mathrm{CVaR}_{0.10}(D) = \mathbb{E}\left[\, D(\theta) \;|\; D(\theta) \leq q_{0.10}(D) \,\right], \qquad \theta \sim p(\theta \,|\, Y), \qquad x^{*} = \arg\max_{x}\; \mathrm{CVaR}_{0.10}\!\left(D(x;\, \theta)\right) \tag{24}
```



```math
P\!\left(\, \tilde{y}_i(x^{*}) \in S_i \;|\; Y \,\right) \geq 0.95 \qquad \mathrm{for\ every\ hard\!-\!limit\ response\ } i \tag{25}
```


- P(osmolality within limit) ≥ 0.95, P(pCO₂ within limit) ≥ 0.95, P(viability at harvest ≥ limit) ≥ 0.95. Thresholds are proposals for team review and are frozen at Round 0.

On a flat ridge, this rule automatically selects the point that performs acceptably even when the model is wrong, rather than the fragile peak of a point-estimate surface. The expected cost is a slightly lower expected titer than the point-estimate optimum; the expected gain is a higher first-pass confirmation rate and a wider usable operating region. That trade is intentional under a robustness-first brief.

## PAR from the posterior upper bound; NOR from measured control noise

- **Conservative where uncertain.** PAR: computed from the local J at the chosen setpoint as in the original — the ellipsoid δᵀ V Σ² Vᵀ δ ≤ 1, half-width 1/σₖ along vₖ — except that σₖ is replaced by its posterior upper credible bound. Where the model is least certain, the proven acceptable range shrinks automatically.

- **Measured, not assumed.** NOR: the factor κ = 0.5 in the original is replaced by a calculation. The empirical distribution of (actual − setpoint) for each factor is estimated from the Round 1 runs (the data model already records both); the NOR is the largest region such that, convolving that measured control-error distribution with the local model, each decision response stays within specification with probability at least 0.95. NOR width therefore reflects measured control capability, not a fixed constant.

In formulas, with Σ_up the diagonal matrix of posterior upper credible bounds σₖ,up from (20) and V̂ the posterior-median recipes at x\*, the PAR used here is the original ellipsoid (the original's rule, ‖Jδ‖₂ ≤ 1) with the uncertain gains replaced by their upper bounds — so the ellipsoid, and each half-width 1/σₖ,up, shrinks exactly where the posterior is widest:


```math
\delta^{\top} \hat{V}\, \Sigma_{\mathrm{up}}^{2}\, \hat{V}^{\top} \delta \;\leq\; 1, \qquad \mathrm{half\!-\!width\ along\ } v_k = \frac{1}{\sigma_{k,\,\mathrm{up}}} \tag{26}
```


And the NOR replaces the original's fixed κ = 0.5 with the measured control error ε = actual − setpoint, whose empirical distribution is estimated from the Round 1 runs. The NOR is the largest region R around x\* such that a run aimed anywhere in R, perturbed by a control error drawn from that empirical distribution and by posterior uncertainty, still meets every specification with probability at least 0.95:


```math
\mathrm{NOR} = \mathrm{largest\ } R:\quad P\!\left(\, \tilde{y}(x^{*} + \delta + \varepsilon) \in S\ \forall\, i \;\;|\;\; Y,\ \varepsilon \sim \mathrm{empirical} \,\right) \geq 0.95, \quad \delta \in R \tag{27}
```


Reporting is unchanged: because operators set factors one at a time, the reported PAR and NOR are the largest upright boxes inside the respective ellipsoids and the design region, with any coupled limit on a stiff mode (for example a joint pH–DO rule) stated explicitly when the box becomes impractically thin.

# 7 Round 2 — Adaptive correction in 3–4 runs, with a quantitative gate

Round 2 is not open-ended augmentation. Before any Round 2 run is placed, the team computes, by posterior simulation, the probability that a confirmation run at the current setpoint passes in both vessels. That single number gates the round:

The gate is a Monte Carlo average over posterior draws: for each draw θ⁽ˢ⁾, simulate a confirmation pair at the current setpoint from the posterior predictive (23) and record whether every decision response meets its criterion in both vessels. The fraction of passing draws is the pass probability:


```math
\hat{p}_{\mathrm{pass}} = \frac{1}{S} \sum_{s=1}^{S} \mathbf{1}\!\left\{\, \tilde{y}_{\mathrm{new}}^{(s)} \ \mathrm{meets\ every\ criterion,\ both\ vessels} \,\right\} \tag{28}
```


| **Simulated pass probability** | **Round 2 action**                                                                        | **Runs**         |
|--------------------------------|-------------------------------------------------------------------------------------------|------------------|
| ≥ 0.80 (proposed threshold)    | Shrink Round 2 or merge into confirmation                                                 | 0–2              |
| 0.50 – 0.80                    | Targeted correction (see below)                                                           | 3–4              |
| \< 0.50                        | Stop and review: model or ranges are wrong; do not spend runs refining a fragile setpoint | 0 (review first) |

When correction runs are placed, they go, in priority order, to: (1) locations where the quadratic and GP models disagree most; (2) the neighbourhood of the candidate setpoint, to sharpen the local J that determines PAR and NOR; and (3) the nearest PAR boundary in the stiff directions. Follow-up is designed in the stiff subspace (coordinates along vₖ), not merely over 2–4 'active factors' — the stiff modes have lower dimension than the active-factor set, so each run carries more decision-relevant information. If the Round 1 bridge flagged a vessel interaction, Round 2 runs are split across vessels instead.

# 8 Round 3 — Confirmation in 2 runs, judged on pre-registered criteria

One run per vessel at the chosen setpoint (2 runs total). Confirmation passes when every decision response falls inside its 95% posterior predictive interval, computed before the runs start. The criteria, the intervals and the setpoint are frozen at the end of Round 2 at the latest; they are not revised after confirmation data arrive.

A failed confirmation does not restart the study. As in the original, it adds runs where the model was wrong — but under this design the failure is also diagnostic: a response outside its predictive interval, with a well-fitting Round 1–2 model, points first to an unmodelled block effect, a control-loop covariate shift, or a vessel effect the bridge under-estimated, in that order of checking.

# 9 Run budget and timeline — head-to-head with the original design

|                    | **Original design**                               | **This alternative**                                                  |
|--------------------|---------------------------------------------------|-----------------------------------------------------------------------|
| Screening          | About 21 runs (DSD, vessel crossed, fake factors) | 16 runs (13-run DSD + 3-run bridge)                                   |
| Follow-up          | Augmentation, open-ended; typically 2–6+ runs     | 3–4 runs, gated; can shrink to 0–2                                    |
| Confirmation       | 2 runs (1 per vessel)                             | 2 runs (1 per vessel)                                                 |
| Total culture runs | About 25–30                                       | About 21 (19–22 depending on the Round 2 gate)                        |
| Culture rounds     | 3 + augmentation rounds                           | 3                                                                     |
| Analysis           | Forward selection; bootstrap CIs                  | Shrinkage posterior; dual quadratic / GP; bootstrap cross-check       |
| Setpoint rule      | Maximize D at the point estimate                  | Maximize CVaR of D, with chance constraints                           |
| PAR / NOR          | PAR from local J; NOR = κ/σₖ, κ = 0.5             | PAR from posterior upper bound of σₖ; NOR from measured control noise |

Timeline follows the original's arithmetic: each block takes about 3 weeks (14-day run plus turnaround), and the total calendar time is governed by how many vessels can run simultaneously — the original's open question 1, which no design can circumvent. If only two or three vessels are available at once, Round 1 itself splits into multiple blocks; the run count is unchanged but the calendar lengthens. Vessel characterization (Phase 0, 1–2 weeks) and design / data-model sign-off (1 week, in parallel) are unchanged, with the Round 0 prior-elicitation session added to the sign-off week.

| *Honest limit: the run saving here is about 15–25%, not a step change. A materially smaller screening design is available only by trusting the priors to place the runs — which this alternative deliberately refuses to do.* |
|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|

# 10 Risks, trade-offs and failure modes

| **Risk / trade-off**                                      | **Effect**                                                                                       | **Mitigation in this design**                                                                                                         |
|-----------------------------------------------------------|--------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------|
| Priors are wrong (flask platform does not transfer)       | Shrinkage fit is pulled toward wrong effects; robust setpoint is conservative in the wrong place | Design placement does not use priors; weakly informative fallback; leave-one-block-out gate before any decision                       |
| Vessel interaction is real but the 3-run bridge misses it | One setpoint is used where two are needed                                                        | Bridge points chosen at diagonal extremes; uncertain bridge result triggers vessel-split Round 2; confirmation is per vessel          |
| Analysis discipline is not maintained                     | Pre-registration and shrinkage degrade into ad-hoc selection; the run saving is spent on rework  | Thresholds frozen at Round 0; pipeline (Section 11) encodes the rules, not analyst discretion                                         |
| Robust setpoint is conservative                           | Expected titer slightly below the point-estimate optimum                                         | Accepted by design under a robustness-first brief; report both setpoints and their posterior D distributions for the decision meeting |
| GP in the stiff subspace overfits at very small n         | False model disagreement misdirects Round 2 runs                                                 | GP is a cross-check only; disagreement must exceed a pre-set tolerance in tolerance units before it directs a run                     |
| Few vessels available simultaneously                      | Round 1 splits into many blocks; calendar lengthens; block variance grows                        | Centre point in every block; block as a random effect; no design removes this constraint                                              |
| Output scaling sᵢ is a judgment (as in the original)      | Stiff / sloppy ranking depends on it                                                             | Agree sᵢ at Round 0; report sensitivity of the stiff set and setpoint to plausible alternative sᵢ                                     |

# 11 What is kept, deliverables, and decisions needed before Round 1

## Kept from the original, unchanged in substance

- The multivariable framing (scaled x → scaled y), SVD of J with stiff / sloppy modes, and canonical analysis of the quadratic surface B.

- The factor set, ranges (subject to vessel characterization), fixed factors, and the response definitions, goals and sᵢ framework.

- The JSON data model with setpoint and actual-value records, daily time series retained, and the reusable Python pipeline — extended, not replaced (see deliverables).

- Blocking with a centre point per block, modelling on actual values, and gas / base totals as covariates.

## Deliverables

- Round 0 prior register (elicited effects, noise ranges, frozen thresholds) and the design tables for Rounds 1–3.

- Python pipeline, extended with: shrinkage (regularized horseshoe) fitting with hierarchical block effects; GP model in the stiff subspace with ARD; posterior simulation for CVaR desirability, chance constraints, pass-probability gating, and posterior PAR / NOR calculation.

- Mode report (vₖ, σₖ with posterior intervals, uₖ), canonical-analysis report, robust setpoint with its posterior D distribution, and PAR / NOR per vessel (or a justified single recipe).

- Confirmation report against the pre-registered predictive intervals, and a synthetic-data test of the full pipeline (known modes) before real runs arrive — as proposed in the original.

## Decisions needed before Round 1

- How many 7 L and 15 L vessels can run simultaneously? (Governs blocking and calendar; unchanged open question.)

- Which vessel is the 'main' vessel for the 13-run DSD, and are the bridge diagonal extremes feasible in the second vessel?

- Platform centre values, response goals, hard limits and sᵢ — plus the Round 0 prior register and the frozen thresholds in Sections 6–8 (CVaR level, chance-constraint levels, pass-probability gate).

- Is product quality part of the optimization objective, or only a constraint? (Unchanged open question; it changes the desirability weights and therefore the robust setpoint.)

**Bottom line for the review meeting.** If the team's binding constraint is calendar rounds and first-pass success, this alternative is the safer buy: fewer runs, no open-ended augmentation, and a setpoint chosen to survive model error. If the binding constraint is maximum expected titer and the team is willing to spend augmentation runs chasing it, the original's point-estimate rule will, on average, aim slightly higher — and fail confirmation more often. The two documents share enough structure that the team can also adopt this design's Sections 5–7 (analysis and decision rules) on top of the original's 21-run screening design, at no extra run cost.

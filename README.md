# Pre-Injury Inhibitory Tone, Not the Acute Injury Response, Predicts Post-Traumatic Epileptogenesis in Larval Zebrafish

**A startle-habituation biomarker tested against its own falsifying alternative, validated by immediate early gene expression and pharmacological seizure threshold**

**ISEF category:** Translational Medical Science (TMED) · *alternates:* Computational Biology & Bioinformatics, Animal Sciences

---

## Hypothesis

> Startle habituation in larval zebrafish is produced by feedforward inhibition onto the Mauthner
> cell, so its decay constant τ indexes inhibitory gain in a defined sensorimotor circuit. Because
> post-traumatic epileptogenesis begins as a loss of inhibitory control, **τ measured around a blast
> injury and interpreted alongside blast dose will identify which individual larvae go on to develop
> spontaneous seizure-like activity.**

The design tests this against the explanations that would make it uninteresting: general sickness,
injury severity alone, and pre-existing trait alone. Each is fitted as its own model under identical
cross-validation, so the hypothesis can fail in a specific and readable way.

**It partly did.** See section 5.

---

## Result in one paragraph

Conversion to spontaneous burst activity at 6.8 dpf is predictable from habituation kinetics at a
cross-validated **AUC of 0.787** and **total out-of-fold accuracy of 70.7%**, with the predictor set
chosen inside the cross-validation and a permutation *p* = 0.001 against a null that reruns the whole
selection procedure. But the predictive signal sits in the **pre-injury** habituation constant plus
blast dose, not in the acute change in habituation that the study was built around. The acute Δτ terms
do not separate converters from non-converters within either dose group, and adding them costs 0.051
AUC. Two assays sharing no instrumentation with the behavioural pipeline confirm that the outcome
itself is a real hyperexcitable state: converter pools carry 73.5% more c-fos transcript than matched
non-converter pools, and larvae that seize under pentylenetetrazol are 17 times more likely to have
converted.

Full numbers with every test statistic, effect size and confidence interval: **[RESULTS.md](RESULTS.md)**

---

## 1. The problem

Post-traumatic epilepsy is the archetypal acquired epilepsy: a person who was not epileptic sustains
a brain injury and, months to years later, begins having spontaneous recurrent seizures. Risk scales
steeply with injury severity.

Between the injury and the first spontaneous seizure sits the **latent period**. During it the network
is being remodelled, inhibitory interneurons are lost or silenced, excitatory circuits sprout, and the
excitation/inhibition set point drifts. By the measures we routinely apply, the patient looks
seizure-free.

That period is the therapeutic opportunity and the scientific bottleneck at once. It is the only
window in which an anti-epileptogenic treatment could work, because it is the window before the
epileptic network exists. But we cannot tell who is in it. Most injured patients never develop
epilepsy, so without a way to enrich a trial for the minority who will, a prevention trial needs an
impractical number of participants and an impractical follow-up. **No validated early biomarker of
epileptogenesis currently exists.**

## 2. Neurobiological rationale

This section is the argument that the dependent variable measures something real. Every design choice
downstream follows from it.

### 2.1 Startle habituation is a readout of inhibition, not of fatigue

Larval zebrafish respond to an abrupt acoustic-vibrational stimulus with a **C-start** escape, a
stereotyped sub-10-millisecond body bend commanded by the **Mauthner cell**, a single identified
reticulospinal neuron in each hindbrain hemisphere. This is one of the best-characterised sensorimotor
circuits in vertebrate neuroscience, which is precisely why it was chosen.

Under repeated stimulation the response habituates. The critical mechanistic point is what habituation
is *not*: it is not depletion or fatigue of the Mauthner cell. Calcium imaging in intact larvae showed
that habituation is driven by **progressively decreased excitability of the M-cell lateral dendrite**,
produced by feedforward inhibitory drive onto that dendrite (Marsden & Granato, 2015). Habituation is
an actively maintained inhibitory process, and it is pharmacologically dissociable from the startle
response itself (Burgess & Granato, 2007; Wolman et al., 2011).

The decay constant τ, trials to habituate, is therefore an index of how rapidly inhibitory gain
accumulates.

> **A longer τ means inhibition builds more slowly: reduced inhibitory gain, a network shifted toward
> excitation.**

### 2.2 Why that is the right variable for epileptogenesis

Post-traumatic epileptogenesis is, at circuit level, a **collapse of the excitation/inhibition
balance**: selective vulnerability and functional silencing of GABAergic interneurons after injury
(Sloviter, 1991), loss of inhibitory control over glutamatergic networks, and aberrant excitatory
reorganisation during the latent period (Hunt, Boychuk & Smith, 2013).

τ and epileptogenesis are not two loosely correlated phenomena. They are two readouts of **the same
underlying quantity**, inhibitory tone. That is what makes τ a mechanistically motivated candidate
rather than a fishing expedition, and it is why the model's coefficients can be sign-checked against
biology instead of merely reported.

### 2.3 Why blast dose must be in the model: the biphasic response

| Dose | τ moves | Mechanism | Interpretation |
|---|---|---|---|
| **Low impact** | **increases** (+3.29 trials at 0.5 h) | Sublethal injury compromises the feedforward inhibition accumulating onto the M-cell dendrite | **Habituation deficit**, disinhibition, the classic hyperexcitable phenotype |
| **High impact** | **decreases** (−1.94 trials at 0.5 h) | Greater energy deposition also depresses the excitatory limb: acute metabolic crisis and depolarisation reduce the startle response itself | **Fatigue, not learning**: the decay is fast because the response never had far to fall |

Both are pathological. Both are evidence of injury. They move the same scalar in opposite directions
(Welch *t* = 9.01, *p* < 0.0001, Cohen's *d* = 1.97).

A model handed Δτ without dose is asked to treat +3 trials and −2 trials as opposite kinds of evidence
when they are the same kind of evidence about two different lesions. Encoding dose resolves the
ambiguity. **This group-level effect is robust and replicated across both datasets analysed.** Whether
it translates into individual-level prediction is a separate question, and section 5 is where it is
answered.

### 2.4 Why a pre-injury baseline session

Every larva is measured **before** injury (t = −1 h) and serves as its own control. Individual
variation in baseline inhibitory tone is large, and an animal that starts nearer its seizure threshold
has less reserve to lose. Using within-animal change removes that variance from the injury signal;
retaining pre-injury τ as its own predictor tests the complementary two-hit idea of susceptibility
plus insult. As it turns out, that second term is the one that carries the signal.

## 3. Variables

| Role | Variable | Operationalisation |
|---|---|---|
| **Independent** | Blast dose | `sham` / `low_impact` / `high_impact`, delivered at 4 dpf |
| | Time since injury | −1 h (pre-injury baseline), 0.5 h, 24 h |
| **Dependent (primary)** | Conversion | `converted` = spontaneous burst activity at 6.8 dpf (0/1) |
| **Dependent (mediator)** | Habituation decay constant τ | Nonlinear fit of `A·exp(−(k−1)/τ) + C` across 30 trials |
| **Molecular outcome** | c-fos expression | `fosab` fold change vs `rpl13a`, 2^−ΔΔCт |
| **Pharmacological outcome** | Seizure threshold | Proportion seizing and latency under 2.5 mM PTZ |
| **Controlled** | Clutch, well position, trial count, stimulus protocol, PTZ concentration, age at injury and at outcome | Identical protocol across all 12 sessions |
| **Blocking factor** | Clutch (A/B/C) | 3 independent clutches on separate days, held out whole in cross-validation |

**One cohort.** All 132 larvae are `zf_*`: every animal contributes habituation sessions and an
outcome, a subset of 81 contributes to a qPCR pool, and a subset of 54 receives a PTZ challenge. This
differs from the pilot design and has a consequence stated throughout: the molecular and
pharmacological analyses **share animals** with the model cohort, so they validate the outcome label
rather than the model's generalisation. Section 6 is explicit about this.

## 4. Methods

### 4.1 Rebuilding the dependent variable

τ is refit from raw trial-level data rather than taken on trust. For each fish-session,
`distance_mm(k) = A·exp(−(k−1)/τ) + C` is fitted across all 30 trials by **nonlinear least squares**
(`scipy.optimize.curve_fit`, bounded, four starting points, best sum-of-squared-error retained).
390/390 sessions converged, median R² = 0.836. The refit reproduces the supplied `decay_constant` at
*r* = 1.0000.

**Log-linearisation was rejected, and the pipeline demonstrates why rather than asserting it.**
Regressing log(y − C) on trial number requires y > C. Once responses reach the habituated floor they
scatter below it: **26.2% of trials are unusable**, and discarding them tilts the fitted slope so far
that correlation with the reference τ collapses from *r* = 1.000 to *r* = −0.252. It inverts. This
diagnostic runs on every execution.

### 4.2 Prediction model

Injured larvae only. One row per fish. **n = 82, 34 converters (41.5%).**

| Predictor | Definition |
|---|---|
| `dose` | high_impact = 1, low_impact = 0 |
| `pre_tau` | τ at t = −1 h |
| `z_dtau_0.5` | τ@0.5 h − τ@−1 h, z-scored **within dose group** |
| `z_dtau_24` | τ@24 h − τ@−1 h, z-scored **within dose group** |

**Model: L2-penalised (ridge) logistic regression in a `StandardScaler` pipeline.** No random forest,
no gradient boosting, no ensemble. At this sample size an ensemble has enough capacity to memorise the
data, its fold-to-fold variance exceeds the effect being measured, and its output cannot be
sign-checked against section 2.

Events per variable is **8.5**, below the conventional floor of 9 to 10 (Peduzzi et al., 1996). The
predictor set is held at four rather than trimmed to satisfy the rule, because each term answers a
distinct pre-specified question. The cost is wider coefficient intervals, and they are reported as
such. The two-predictor model the data actually support clears the rule comfortably at EPV 17.

**Validation, three separate fixes for three separate problems:**

1. **`GroupKFold` on clutch for the outer split** (leave-one-clutch-out), never random. Clutches are
   sibling groups run on separate days. Quantified in section 5.4: a random split inflates AUC by
   0.091.
2. **Nested cross-validation for the penalty strength C**, selected inside each outer training set
   only, never on the folds reported (Varma & Simon, 2006).
3. **Model selection inside the cross-validation.** Quoting the best row of an ablation table is
   selection on the test folds. The choice of predictor set is instead treated as a hyper-parameter
   tuned on inner folds, so the reported AUC is honest about the selection.

**Permutation test:** labels shuffled **within clutch**, preserving each clutch's conversion rate, and
the entire nested procedure rerun 1,000 times.

## 5. Results

### 5.1 Groups start equal

At the pre-injury baseline, τ does not differ between groups (ANOVA *F*(2,126) = 2.06, *p* = 0.132,
η² = 0.032; Kruskal–Wallis *p* = 0.153), nor does baseline locomotion (*p* = 0.087). Randomisation
held.

### 5.2 The ablation: where the hypothesis breaks

| Model | Question it answers | Mean fold AUC | SD | Fold range |
|---|---|---|---|---|
| (a) `baseline_locomotion` only | Is it just sickness? | 0.392 | 0.066 | 0.322–0.453 |
| (b) `dose` only | Is it just injury severity? | 0.623 | 0.044 | 0.592–0.674 |
| **(c) `dose + pre_tau`** | Is it a pre-existing trait? | **0.788** | 0.061 | 0.722–0.841 |
| (d) `dose + z_dtau` | Is it the injury response? | 0.627 | 0.097 | 0.547–0.735 |
| (e) full 4-predictor model | *(pre-specified)* | 0.738 | 0.083 | 0.667–0.829 |

- **Not sickness.** Baseline locomotion alone is below chance out of fold.
- **Not severity alone.** Dose alone reaches 0.623, real but modest.
- **The acute injury response does not carry individual-level signal.** Model (d) barely improves on
  dose alone, and Δτ fails to separate converters from non-converters within either dose group
  (low impact *p* = 0.475, high impact *p* = 0.631).
- **Adding the Δτ terms to `pre_tau` costs 0.051 AUC.** Model (e), the pre-specified model, is beaten
  by the simpler model (c).

**The central hypothesis is not supported in the form it was posed.** What predicts conversion is the
pre-injury habituation constant together with dose. Reporting this is the point of running an ablation
rather than fitting one model and describing it.

The group-level Δτ effect is not in doubt (section 2.3, *d* = 1.97). What fails is the step from a
group difference to an individual prediction: within a dose group, the size of a larva's acute Δτ does
not tell you whether that larva will convert.

### 5.3 Honest performance

Because model selection is itself a source of optimism, the predictor set is chosen on inner folds
only:

| Held-out clutch | Selected by inner CV | Outer-fold AUC |
|---|---|---|
| clutch_A | `e_full` | 0.717 |
| clutch_B | `c_dose_pretau` | 0.802 |
| clutch_C | `c_dose_pretau` | 0.841 |

| Metric | Value |
|---|---|
| **Mean fold AUC (honestly selected)** | **0.787** (SD 0.064) |
| Pooled out-of-fold AUC | 0.746 [0.629, 0.848] |
| **Total out-of-fold accuracy** | **70.7%** |
| Balanced accuracy | 69.0% |
| Sensitivity / Specificity | 58.8% / 79.2% |
| Permutation *p* (whole selection procedure) | **0.001** (floor at 1,000 permutations; null mean 0.487) |
| Optimism avoided by nesting the selection | +0.002 AUC |

### 5.4 Leakage quantification

| Split | Mean fold AUC |
|---|---|
| 5-fold **random** (leaky, for comparison only) | 0.828 |
| Leave-one-clutch-out (**the honest estimate**) | 0.738 |

Ignoring clutch inflates AUC by **0.091**, five times larger than in the pilot cohort. Any figure from
a random split on this design should be discounted accordingly.

## 6. Orthogonal validation

*Two assays, neither sharing instrumentation with the behavioural pipeline, test whether `converted`
denotes a real hyperexcitable state. Both are run on larvae drawn from the same cohort as the model,
so they validate the **outcome label**, not the model's generalisation to unseen animals. That
distinction is maintained throughout.*

### 6.1 Why it is needed

Everything in section 5 derives from one number per trial: how far a larva swam. If that pipeline
contains a systematic artefact, **no amount of internal cross-validation would reveal it**, because
every fold inherits the same artefact. Held-out clutches protect against overfitting; they cannot
protect against a measurement that is measuring the wrong thing. The only way to test that is to leave
the modality.

### 6.2 c-fos: the molecular axis

`fosab` is the zebrafish orthologue of *c-fos*, the canonical immediate early gene. Sustained
depolarisation raises intracellular Ca²⁺, driving CaMK- and MAPK/ERK-dependent phosphorylation of CREB
and transcription from the *fos* promoter within roughly 15 to 30 minutes (Sheng & Greenberg, 1990).
c-fos transcript is a molecular integrator of recent network activity and the standard readout for
mapping seizure-recruited circuits, including in the original characterisation of chemically induced
seizures in larval zebrafish (Baraban et al., 2005). Normalisation is against `rpl13a`, a validated
zebrafish reference gene (Tang et al., 2007), by the 2^−ΔΔCт method (Livak & Schmittgen, 2001).

**The pools are matched, and the matching is unbalanced by biology.** There are 27 pools of 3 larvae,
labelled by realised outcome. Because conversion is rare in sham, sham cells contain **no converter
pool at all**:

| Group | Converter pools per clutch | Non-converter pools per clutch |
|---|---|---|
| sham | 0 | 2 |
| low_impact | 1 | 2 |
| high_impact | 2 | 2 |

A balanced nine-pair design is therefore impossible. The pairing the data support is **6 matched
cells**, the injured group × clutch combinations, with replicate pools averaged within a cell before
pairing. Treating the 27 pools as independent would inflate the degrees of freedom and ignore the
matching.

| Contrast | Pairs | Mean Δ (converter − non-converter) | 95% CI | *t* | *p* | Cohen's *dz* | Wilcoxon *p* |
|---|---|---|---|---|---|---|---|
| Fold-change scale | 6 | +0.955 | [+0.296, +1.614] | *t*(5) = 3.73 | **0.014** | 1.52 | 0.031 |
| log₂ scale | 6 | +0.795 | [+0.307, +1.283] | *t*(5) = 4.19 | **0.009** | 1.71 | 0.031 |

Converter pools carry **73.5% more c-fos transcript** relative to `rpl13a` (geometric mean ratio 1.735,
95% CI [1.237, 2.434], back-transformed from the paired log₂ analysis). Direction is consistent in **6 of 6** pairs (sign test
*p* = 0.031).

**Specificity control.** Since sham cannot be paired, the control question is turned around: if injury
alone raised c-fos, injured **non-converter** pools would sit above sham non-converter pools. They do
not (1.256 vs 1.181, Welch *t* = 0.56, *p* = 0.607, *d* = 0.39). **A null result here is the desired
one.** The c-fos elevation tracks which larvae converted, not which larvae were hit. That is what makes
the assay informative about epileptogenesis rather than about injury exposure.

### 6.3 PTZ: the pharmacological axis

Pentylenetetrazol is a non-competitive GABAₐ receptor antagonist that binds at the picrotoxin site,
reduces chloride conductance and removes inhibitory brake from the network (Baraban et al., 2005). It
is the standard probe of **seizure threshold**.

| Group | Seized / n | Proportion | Wilson 95% CI | Median latency (s) |
|---|---|---|---|---|
| Sham | 5/17 | 0.29 | [0.13, 0.53] | 1800 |
| Low impact | 11/15 | 0.73 | [0.48, 0.89] | 902 |
| High impact | 19/22 | 0.86 | [0.67, 0.95] | 589 |

χ²(2) = 14.30, *p* = 0.0008, Cramér's V = 0.515, n = 54. Post-hoc power for the sham-versus-injured
contrast is **96%**, so unlike the pilot cohort this probe is adequately powered. It remains a
secondary outcome for a different reason: it is a group-level comparison, and the primary claim of this
project is about individual animals. The low-versus-high contrast is still not powered (*p* = 0.408).

**Is the challenge a confound?** PTZ is a proconvulsant given to a subset of the same larvae that
supply the outcome, so the obvious worry is that the drug caused the outcome. It did not. Conversion
rates are effectively identical between challenged and unchallenged larvae in every group (sham
*p* = 1.000, low impact *p* = 0.737, high impact *p* = 1.000).

**Seizure threshold and conversion coincide at the level of the individual animal.** Larvae that seized
under PTZ were far more likely to have converted (Fisher exact **OR = 17.0**, *p* = 0.0019), and
converters reached their first seizure much sooner (median 405 s vs 1775 s, Mann–Whitney *p* = 0.0008).
This was not pre-specified and is exploratory, but it is the strongest evidence in the project that
`converted` denotes a real hyperexcitable state rather than a scoring artefact. It says nothing about
whether that state is *predictable in advance*, which is section 5's job.

## 7. Conclusion

Conversion is predictable at **AUC 0.787** and **70.7% total accuracy** under leave-one-clutch-out
nested cross-validation with the model chosen inside the folds, permutation *p* = 0.001. Sickness and
injury severity were each tested as alternative explanations and each falls short.

The signal lies in **pre-injury inhibitory tone plus blast dose**, not in the acute injury response the
study was designed to test. Blast does perturb the circuit strongly and in dose-dependent opposite
directions, but that group-level effect does not carry individual-level predictive information here.

This matters for the framing. A variable measured before an injury is a **susceptibility marker**, not
an acute readout, and it could not be used to triage a patient after an injury that has already
happened. The clinical use case the project began from is not supported by this cohort. What is
supported is narrower and still worth having: **individual variation in inhibitory gain, measurable
non-invasively and cheaply, predicts which animals undergo epileptogenesis after a standardised
insult**, and the epileptogenic state itself is confirmed by two independent modalities.

The behavioural assay costs **$0.09 and 0.9 operator-minutes per animal**.

## 8. Limitations

- **The central hypothesis is not supported in the form posed.** The acute Δτ terms add nothing and
  slightly degrade prediction. This is reported as the primary finding it is, not buried.
- **n = 82 with 34 events, EPV = 8.5**, below the conventional floor. Intervals are wide.
- **Three clutches means three outer folds.** Fold spread is estimated from three numbers.
- **6 injured larvae were dropped** for a missing session. Complete-case analysis; attrition is not
  modelled.
- **Conversion is a behavioural proxy.** Not electrographically confirmed epilepsy, though the PTZ
  concordance in 6.3 supports it. Tectal field recording is the confirmatory measurement.
- **The c-fos and PTZ analyses share animals with the model cohort.** They validate the outcome label,
  not generalisation.
- **The c-fos assay is bulk and pooled.** 6 pairs, 3 larvae per pool. It cannot localise the signal to
  a cell type, and cannot separate loss of parvalbumin-positive interneuron function from increased
  glutamatergic drive.
- **Cell-type resolution is absent.** The inhibitory-failure account in section 2 is inference from
  circuit behaviour, not direct observation in these animals.
- **One species, one age, one injury model.**

## 9. Future work

1. **Electrophysiological confirmation.** Tectal field recordings to establish that conversion is
   electrographic seizure activity.
2. **Test the susceptibility finding directly.** If `pre_tau` marks inhibitory reserve, larvae in the
   top and bottom baseline-τ quartiles should differ in conversion rate under identical injury. That
   is a clean, pre-registrable experiment this dataset motivates.
3. **Denser acute sampling.** Δτ was sampled at 0.5 h and 24 h. The pilot cohort had 1 h and 5 h as
   well and still showed no individual-level signal, but a finer trajectory would settle it.
4. **Cell-type resolution.** Transgenic GABAergic reporter lines to test the interneuron-failure
   account directly.
5. **Pharmacological rescue.** If τ indexes inhibitory gain, an agent restoring inhibition during the
   latent period should reduce conversion and renormalise τ. That is the experiment that moves this
   from biomarker to mechanism.
6. **External replication.** Three clutches is the binding constraint on every uncertainty estimate
   here.

## 10. ISEF compliance and ethics

- **Vertebrate animal research.** Under ISEF rules zebrafish are not considered vertebrate animals
  until 72 hours post-fertilisation. This study injures larvae at 4 dpf and follows them to 6.8 dpf, so
  it **is** regulated vertebrate animal research. It requires prior SRC/IACUC approval, Form 5B, a
  Qualified Scientist, and a Regulated Research Institution. **Verify current-year rules with your
  fair's SRC before starting; do not rely on this summary.**
- **Humane endpoints.** Blast injury and PTZ challenge both cause more than momentary distress and must
  be covered explicitly in the approved protocol, with anaesthesia and euthanasia methods specified.
- **Reduction.** 132 larvae total, 12 recording sessions, 4 lost. Terminal assays are drawn from the
  followed cohort rather than requiring separate animals, which reduces total use.
- **Data integrity.** Every reported statistic is machine-generated into
  [`results/all_statistics.csv`](results/all_statistics.csv). Nothing in RESULTS.md is hand-transcribed.
  The seed is fixed and reported; reruns are byte-identical.

## 11. Running it

```bash
pip install -r requirements.txt
python run_all.py
```

Options:

```bash
python run_all.py --permutations 200   # faster run (default 1000)
python run_all.py --skip-permutation   # skip the permutation tests
python run_all.py --seed 123           # override the random seed
```

The full run takes a few minutes, most of it the two 1,000-iteration permutation tests. Each iteration
reruns the entire nested cross-validation, and for the selection test, the model selection too.

### Data

`tbidataset.xlsx` at the repository root, six sheets: `habituation_trials`, `fish_features`,
`outcomes`, `cfos_pools`, `ptz_challenge`, `session_log`. To point the pipeline elsewhere, change
`DATA_XLSX` in [`src/config.py`](src/config.py). The loader aliases several column names, so
`delta_delta_ct`/`delta_ddct` and `risk_pool`/`pool_type` are both accepted.

### Repository layout

```
run_all.py                  runner, executes every step in order
src/
  config.py                 seed, paths, design constants, plot style
  statsbook.py              append-only ledger; everything reported passes through it
  io_data.py                workbook loading, sheet validation, column aliasing
  habituation.py            STEP 1, nonlinear curve fitting and agreement checks
  features.py               design matrix, model specifications
  modeling.py               ridge logistic, nested CV, nested model selection, permutation
  step2_prediction.py       STEP 2 driver, ablation, leakage, coefficients
  step3_cfos.py             STEP 3, paired c-fos pools
  step4_descriptive.py      STEP 4, descriptives, operations, PTZ
  plotting.py               all figures, 300 dpi PNG
  report.py                 assembles RESULTS.md from the ledger
results/
  figures/*.png             13 figures at 300 dpi
  tables/*.csv              intermediate tables
  all_statistics.csv        every reported statistic, one row each
RESULTS.md                  the written report
```

### Figures

| | |
|---|---|
| `fig01` | Example nonlinear fits with fitted offset |
| `fig02` | τ agreement: refit vs supplied, Bland–Altman, log-linearisation failure |
| `fig03` | **ROC curves for all five ablation models on one axes** |
| `fig04` | **Permutation null histogram with observed AUC marked** |
| `fig05` | Confusion matrix, calibration curve, score separation |
| `fig06` | Coefficients with clutch-clustered bootstrap CIs |
| `fig07` | Per-fold AUC by model |
| `fig08` | **6 matched c-fos pairs, plus the specificity control** |
| `fig09` | Habituation curves by group and timepoint, including t = −1 |
| `fig10` | τ by timepoint and group with error bars |
| `fig11` | Converter vs non-converter τ trajectories, split by dose |
| `fig12` | Operational metrics: throughput, operator time, cost, attrition |
| `fig13` | PTZ proportions with Wilson CIs, and latency distributions |

### Reproducibility

Seed `20260809`, set in [`src/config.py`](src/config.py) and reported at the top of every run and in
RESULTS.md. Assumption checks print tagged `[PASS]` or `[FLAG]`, including the ones that legitimately
flag: EPV below the conventional floor, and non-normal baseline τ in the injured groups (backed by
Kruskal–Wallis, which agrees).

---

## Background literature

Reading behind the neurobiological claims. These support the framing; they are not citations of results
generated here.

1. Annegers JF, Hauser WA, Coan SP, Rocca WA (1998). A population-based study of seizures after traumatic brain injuries. *New England Journal of Medicine* 338:20–24.
2. Baraban SC, Taylor MR, Castro PA, Baier H (2005). Pentylenetetrazole induced changes in zebrafish behavior, neural activity and c-fos expression. *Neuroscience* 131:759–768.
3. Burgess HA, Granato M (2007). Sensorimotor gating in larval zebrafish. *Journal of Neuroscience* 27:4984–4994.
4. Hunt RF, Boychuk JA, Smith BN (2013). Neural circuit mechanisms of post-traumatic epilepsy. *Frontiers in Cellular Neuroscience* 7:89.
5. Livak KJ, Schmittgen TD (2001). Analysis of relative gene expression data using real-time quantitative PCR and the 2^−ΔΔCт method. *Methods* 25:402–408.
6. Marsden KC, Granato M (2015). In vivo Ca²⁺ imaging reveals that decreased dendritic excitability drives startle habituation. *Cell Reports* 13:1733–1740.
7. Peduzzi P, Concato J, Kemper E, Holford TR, Feinstein AR (1996). A simulation study of the number of events per variable in logistic regression analysis. *Journal of Clinical Epidemiology* 49:1373–1379.
8. Sheng M, Greenberg ME (1990). The regulation and function of c-fos and other immediate early genes in the nervous system. *Neuron* 4:477–485.
9. Sloviter RS (1991). Permanently altered hippocampal structure, excitability, and inhibition after experimental status epilepticus in the rat: the 'dormant basket cell' hypothesis. *Hippocampus* 1:41–66.
10. Tang R, Dodd A, Lai D, McNabb WC, Love DR (2007). Validation of zebrafish (*Danio rerio*) reference genes for quantitative real-time RT-PCR normalization. *Acta Biochimica et Biophysica Sinica* 39:384–390.
11. Varma S, Simon R (2006). Bias in error estimation when using cross-validation for model selection. *BMC Bioinformatics* 7:91.
12. Wolman MA, Jain RA, Liss L, Granato M (2011). Chemical modulation of memory formation in larval zebrafish. *PNAS* 108:15468–15473.

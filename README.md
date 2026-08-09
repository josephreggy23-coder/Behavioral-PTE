# Predicting Post-Traumatic Epileptogenesis from Acute Startle-Habituation Kinetics in a Larval Zebrafish Blast Injury Model

**A behavioural biomarker of inhibitory circuit function, validated against immediate early gene expression**

**ISEF category:** Translational Medical Science (TMED) · *alternates:* Computational Biology & Bioinformatics, Animal Sciences

> **Headline result.** A four-predictor ridge logistic model built on how quickly the acoustic
> startle response habituates — measured before injury and again at 30 minutes and 24 hours after
> blast — separates injured larvae that later develop spontaneous seizure-like burst activity from
> those that do not. **Cross-validated AUC 0.833** (95% CI 0.736–0.916), **total out-of-fold accuracy
> 76.5%**, permutation *p* = 0.001 against a null that reruns the entire nested cross-validation
> 1,000 times. An independent qPCR assay on a separate cohort of larvae agrees.

Full numbers, with every test statistic, effect size and confidence interval: **[RESULTS.md](RESULTS.md)**

---

## 1. The problem

Post-traumatic epilepsy is the archetypal *acquired* epilepsy: a person who was not epileptic
sustains a brain injury and, months to years later, begins having spontaneous recurrent seizures.
It accounts for a substantial share of epilepsy in young adults, and risk scales steeply with injury
severity — Annegers and colleagues' population study put the 5-year risk after severe head injury an
order of magnitude above the population baseline, while mild injury carried a much smaller excess.

Between the injury and the first spontaneous seizure sits the **latent period**. During it the
network is being remodelled — inhibitory interneurons are lost or silenced, excitatory circuits
sprout and rewire, and the excitation/inhibition set point drifts — but the animal or patient looks,
by the only measures we routinely apply, seizure-free.

That latent period is the therapeutic opportunity and the scientific bottleneck at once:

- **The opportunity.** It is the only window in which an anti-epileptogenic treatment could plausibly
  work, because it is the window before the epileptic network exists.
- **The bottleneck.** We cannot tell, at the time of injury, who is in that window. Most injured
  patients never develop epilepsy. Without a way to enrich a trial for the minority who will,
  a prevention trial needs an impractical number of participants and an impractical follow-up.

**No validated early biomarker of epileptogenesis currently exists.** That is the gap this project
addresses.

## 2. Research question

> Does an acute, non-invasive behavioural readout of **inhibitory circuit function**, measured within
> 24 hours of a blast injury, predict which individual animals will go on to develop spontaneous
> seizure-like activity?

## 3. Hypotheses

**H₁ (primary).** The *change* in startle-habituation rate from each animal's own pre-injury
baseline, measured at 30 minutes and 24 hours post-blast, predicts conversion to spontaneous burst
activity at 6 days post-fertilisation, at above-chance accuracy under cross-validation that holds out
entire clutches.

**H₂ (dose as moderator).** Low and high blast dose displace habituation rate in **opposite**
directions. Consequently a model given the change in habituation rate *without* blast dose will
perform near chance, because the two signals cancel.

**H₃ (orthogonal validation).** Larvae binned as high-risk on behaviour alone will show elevated
expression of the immediate early gene *fosab* (c-fos) relative to behaviourally low-risk
pool-mates — the molecular signature of a chronically over-active network.

**H₀.** Conversion is unrelated to habituation kinetics and is explained by general sickness
(reduced baseline locomotion), by injury severity alone, or by a pre-existing trait.

---

## 4. Neurobiological rationale

This section is the argument that the dependent variable measures something real. Every design
choice downstream follows from it.

### 4.1 Why startle habituation is a readout of inhibition, not of fatigue

Larval zebrafish respond to an abrupt acoustic-vibrational stimulus with a **C-start** escape — a
stereotyped, sub-10-millisecond body bend. The short-latency C-start is commanded by the **Mauthner
cell**, a single identified reticulospinal neuron in each hindbrain hemisphere. This is one of the
best-characterised sensorimotor circuits in vertebrate neuroscience, which is the reason to use it.

Under repeated stimulation the response habituates. The critical mechanistic point is **what
habituation is not**: it is not depletion or fatigue of the Mauthner cell. Calcium imaging in intact
larvae showed that habituation is driven by **progressively decreased excitability of the M-cell's
lateral dendrite**, produced by feedforward inhibitory drive onto that dendrite (Marsden & Granato,
2015). Habituation is an actively maintained inhibitory process, and it is pharmacologically
dissociable — habituation rate can be modulated without altering the startle response itself
(Burgess & Granato, 2007; Wolman et al., 2011).

So the decay constant of the habituation curve — τ, "trials to habituate" — is an index of how
rapidly inhibitory gain accumulates in a defined circuit.

> **A longer τ means inhibition builds more slowly: reduced inhibitory gain, a network shifted
> toward excitation.**

### 4.2 Why that is the right variable for epileptogenesis

Post-traumatic epileptogenesis is, at circuit level, a **collapse of the excitation/inhibition
balance**. The mechanisms are well documented in mammalian models: selective vulnerability and
functional silencing of GABAergic interneurons after injury (the "dormant basket cell" account of
inhibitory failure, Sloviter, 1991), loss of inhibitory control over glutamatergic networks, and
aberrant excitatory reorganisation during the latent period (Hunt, Boychuk & Smith, 2013).

τ and epileptogenesis are therefore not two loosely correlated phenomena. They are two readouts of
**the same underlying quantity** — inhibitory tone. That is what makes τ a mechanistically motivated
biomarker candidate rather than a fishing expedition, and it is why the model's coefficients can be
sign-checked against biology instead of merely reported.

### 4.3 Why blast dose must be a predictor — the biphasic response

This is the single most important design point, and it is the one an analyst unfamiliar with the
biology would get wrong.

| Dose | τ moves | Mechanism | Interpretation |
|---|---|---|---|
| **Low impact** | **↑ increases** (+3.74 trials at 0.5 h) | Sublethal injury preferentially compromises the feedforward inhibition accumulating onto the M-cell dendrite | **Habituation deficit** — disinhibition, the classic post-injury hyperexcitable phenotype |
| **High impact** | **↓ decreases** (−2.72 trials at 0.5 h) | Greater energy deposition also depresses the excitatory limb — acute metabolic crisis and depolarisation reduce the startle response itself | **Fatigue, not learning** — the decay is fast because the response never had far to fall |

Both are pathological. Both are evidence of injury. **They move the same scalar in opposite
directions** (Welch *t* = 10.00, *p* < 0.0001, Cohen's *d* = 2.19).

A model handed Δτ without dose is being asked to treat +4 trials and −3 trials as opposite kinds of
evidence, when they are the same kind of evidence about two different lesions. Pooled across dose the
two shifts partially cancel and the signal collapses toward chance. Encoding dose tells the model
which direction is pathological for a given animal. **This is why `dose` is a required predictor**,
and the ablation table in §7 is the empirical demonstration.

Note the coefficient on dose is the one whose confidence interval includes zero. That is expected and
correct: dose enters as a **moderator**, not as a main effect.

### 4.4 Why a pre-injury baseline session

Every fish is measured **before** injury (t = −1 h) and serves as its own control. Individual
variation in baseline inhibitory tone is large — visible in the figure below as the spread of
pre-injury τ within every group — and an animal that starts nearer its seizure threshold has less
reserve to lose. Using within-animal change removes that variance from the injury signal, and
retaining pre-injury τ as its own predictor tests the complementary "two-hit" idea: susceptibility
plus insult. The data support both terms carrying independent weight.

### 4.5 Why larval zebrafish

Conserved GABAergic and glutamatergic pharmacology; established seizure models with electrographic
validation; a single identified command neuron for the behaviour being measured; and sufficient
throughput (133 animals, 650 recording sessions, ~9.7 operator-hours total) to power an
individual-level prediction study, which is not feasible at rodent cost.

---

## 5. Variables

| Role | Variable | Operationalisation |
|---|---|---|
| **Independent** | Blast dose | `sham` / `low_impact` / `high_impact` |
| | Time since injury | −1 h (pre-injury baseline), 0.5, 1, 5, 24 h |
| **Dependent (primary)** | Conversion | `converted` = spontaneous burst activity at 6 dpf (0/1) |
| **Dependent (mediator)** | Habituation decay constant τ | Nonlinear fit of `A·exp(−(k−1)/τ) + C` to distance travelled across 30 trials |
| **Orthogonal outcome** | c-fos expression | `fosab` fold change vs `rpl13a`, 2^−ΔΔCт |
| **Secondary outcome** | Seizure threshold | Proportion seizing under fixed 2.5 mM PTZ |
| **Controlled** | Clutch, well position, trial count, stimulus protocol, PTZ concentration, session structure | Identical protocol across all 15 sessions; clutch modelled explicitly |
| **Blocking factor** | Clutch (A/B/C) | 3 independent clutches on separate days; held out whole in cross-validation |

**Cohorts are kept separate by design.** `fish_*` larvae are followed to outcome and are the only
animals in the prediction model. `cf_*` larvae are sacrificed for qPCR and appear only in the
orthogonal validation. `ptz_*` larvae are used only for the threshold probe. No animal contributes to
more than one arm — so §8 is genuinely independent of §7.

---

## 6. Methods

### 6.1 Rebuilding the dependent variable (Step 1)

τ is refit from raw trial-level data rather than taken on trust. For each fish-session,
`distance_mm(k) = A·exp(−(k−1)/τ) + C` is fitted across all 30 trials by **nonlinear least squares**
(`scipy.optimize.curve_fit`, bounded, four starting points per session, best sum-of-squared-error
retained). 650/650 sessions converged; median R² = 0.796. The refit reproduces the supplied
`decay_constant` at *r* = 1.0000 with a bias of −0.001 trials.

**Log-linearisation was rejected, and the pipeline shows why rather than asserting it.** Subtracting
an estimated offset and regressing log(y − C) on trial number requires y > C. Once responses reach
the habituated floor they scatter *below* the offset: **28.1% of trials are unusable**, and
discarding them tilts the fitted slope so far that **4 of 650 sessions return a negative τ** — the
sign inverts, and a hyperexcitable fish is scored as hypo-excitable. Correlation with the reference
τ falls from *r* = 1.000 to *r* = −0.229. This diagnostic runs on every execution.

### 6.2 Prediction model (Step 2)

Injured fish only; sham dropped. One row per fish. **n = 81, 36 converters (44.4%).**

Exactly four predictors:

| Predictor | Definition |
|---|---|
| `dose` | high_impact = 1, low_impact = 0 |
| `pre_tau` | τ at t = −1 h (pre-injury baseline) |
| `z_dtau_0.5` | τ@0.5 h − τ@−1 h, z-scored **within dose group** |
| `z_dtau_24` | τ@24 h − τ@−1 h, z-scored **within dose group** |

**Model: L2-penalised (ridge) logistic regression in a `StandardScaler` pipeline.** No random forest,
no gradient boosting, no ensemble. At n = 81 with 36 events an ensemble has enough capacity to
memorise the sample, its fold-to-fold variance exceeds the effect being measured, and its output
cannot be sign-checked against the biology in §4. A penalised linear model yields four interpretable
coefficients.

**Four predictors is a ceiling, not a starting point.** 36 events ÷ 4 predictors = **9.0 events per
variable**, at the accepted minimum for logistic regression (Peduzzi et al., 1996). A fifth predictor
would push below it.

**Validation — two fixes for two different problems, applied together:**

1. **`GroupKFold` on clutch for the outer split** (leave-one-clutch-out). Never random. Clutches are
   sibling groups run on separate days; a random split puts siblings on both sides of the partition
   and the model learns clutch identity. Quantified in §7: ignoring clutch inflates AUC by +0.016.
2. **Nested cross-validation for the penalty strength.** C is selected inside each outer training
   set only, never on the folds that are reported. Plain k-fold tuning is optimistically biased at
   this sample size (Varma & Simon, 2006).
3. **Permutation test.** Labels shuffled **within clutch** — preserving each clutch's conversion rate,
   the conservative null — and the *entire* nested cross-validation rerun 1,000 times.
4. **Per-fold AUC and spread reported**, not only the mean.

The within-dose z-scoring is itself a data-dependent transform, so it is implemented as a pipeline
step fitted on training folds only (`modeling.WithinDoseZScorer`) rather than applied to the whole
dataset up front. A whole-dataset version is reported as a sensitivity check.

---

## 7. Results

### 7.1 Groups start equal

At the pre-injury baseline, τ does not differ between groups (one-way ANOVA *F*(2,130) = 0.543,
*p* = 0.583, η² = 0.008; Kruskal–Wallis *p* = 0.770), nor does baseline locomotion (*p* = 0.122).
Randomisation held. Everything that follows is a change from a common starting point.

### 7.2 The nested ablation — the scientific argument

Same model class, same cross-validation scheme, different inputs.

| Model | Question it answers | Mean fold AUC | SD | Fold range |
|---|---|---|---|---|
| (a) `baseline_locomotion` only | Is it just sickness? | 0.416 | 0.090 | 0.323–0.503 |
| (b) `dose` only | Is it just injury severity? | 0.440 | 0.079 | 0.354–0.509 |
| (c) `dose + pre_tau` | Is it a pre-existing trait? | 0.677 | 0.132 | 0.588–0.830 |
| (d) `dose + z_dtau` | Is it the injury response? | 0.621 | 0.043 | 0.575–0.659 |
| **(e) full 4-predictor model** | — | **0.849** | 0.080 | 0.758–0.909 |

- **Not sickness.** Baseline locomotion alone is at or below chance out of fold.
- **Not severity.** Dose alone is at chance across clutches — the dose effect is not stable between
  them, exactly as expected for a moderator rather than a main effect.
- **Not trait alone, not response alone.** Each carries real information; neither reaches the full
  model.
- **The full model gains +0.172 AUC over the better of (c) and (d)**, and beats every ablation in
  every one of the three clutch folds. Trait and response carry **complementary** information.

### 7.3 Performance of the full model

| Metric | Value |
|---|---|
| Pooled out-of-fold AUC | **0.833** (95% CI 0.736–0.916) |
| Mean fold AUC | 0.849 (SD 0.080; folds 0.758 / 0.909 / 0.880) |
| **Total accuracy (out-of-fold, threshold 0.50)** | **76.5%** — 62 of 81 fish correct |
| Balanced accuracy | 75.6% |
| Sensitivity / Specificity | 66.7% / 84.4% |
| PPV / NPV | 77.4% / 76.0% |
| Brier score | 0.173 (0.25 = uninformative) |
| Permutation *p* | **0.001** (1,000 reruns of the full nested CV; null mean 0.469) |

Confusion matrix, out-of-fold:

|  | Predicted non-converter | Predicted converter |
|---|---|---|
| **Non-converter** | 38 | 7 |
| **Converter** | 12 | 24 |

The model is more conservative than it is sensitive at this threshold — it misses 12 converters and
false-alarms on 7. For enriching a prevention trial that trade-off is the right way round, but it is
a real limitation, and it is stated as one.

### 7.4 Leakage quantification

| Split | Mean fold AUC |
|---|---|
| 5-fold **random** (leaky — reported only for comparison) | 0.865 |
| Leave-one-clutch-out (**the honest estimate**) | 0.849 |

Ignoring clutch inflates AUC by **+0.016**. Any figure from a random split on this design should be
discounted accordingly.

### 7.5 Coefficients

| Predictor | β (per SD) | 95% CI | Neurobiological reading |
|---|---|---|---|
| `pre_tau` | +0.331 \* | [+0.200, +0.430] | Predisposition — lower baseline inhibitory reserve |
| `z_dtau_24` | +0.257 \* | [+0.133, +0.362] | Failure to renormalise by 24 h → latent period |
| `z_dtau_0.5` | +0.168 \* | [+0.029, +0.290] | Acute disruption of inhibitory gain |
| `dose` | +0.075 | [−0.081, +0.228] | Moderator, not main effect — CI includes zero, as expected |

\* clutch-clustered bootstrap CI excludes zero (2,000 resamples).

---

## 8. Orthogonal validation

*This section is deliberately separate from §7. It is not a robustness check on the model — it is an
independent test of the same hypothesis, using different animals and a different measurement
modality.*

### 8.1 Why it is needed

Everything in §7 derives from one number per trial: how far a larva swam. If that pipeline contains
a systematic artefact — a tracking bias, a plate-position effect, a curve-fitting quirk — **no amount
of internal cross-validation would reveal it**, because every fold inherits the same artefact. Held-out
clutches protect against overfitting; they cannot protect against a measurement that is measuring the
wrong thing.

The only way to test that is to leave the modality entirely.

### 8.2 What is measured

The `cf_*` cohort was sacrificed for quantitative PCR and **never contributed a row to the prediction
model**. Larvae were binned high-risk or low-risk on behaviour and pooled 4 per pool.

**`fosab` is the zebrafish orthologue of *c-fos*, the canonical immediate early gene.** Sustained
neuronal depolarisation raises intracellular Ca²⁺, which drives CaMK- and MAPK/ERK-dependent
phosphorylation of CREB and transcription from the *fos* promoter within roughly 15–30 minutes
(Sheng & Greenberg, 1990). c-fos transcript is therefore a molecular integrator of recent network
activity, and it is the standard readout for mapping seizure-recruited circuits — including in the
original characterisation of chemically induced seizures in larval zebrafish (Baraban et al., 2005).
Normalisation is against `rpl13a`, validated as a stable zebrafish reference gene (Tang et al.,
2007), by the 2^−ΔΔCт method (Livak & Schmittgen, 2001).

**The prediction is specific and falsifiable:** if the behavioural risk score tracks genuine network
hyperexcitability rather than an artefact of swim tracking, high-risk pools should carry more c-fos
transcript than their low-risk pool-mates. Behaviour and transcription share no instrumentation, no
analyst, and no fish.

### 8.3 Statistical treatment — the pools are paired

The 18 pools are **not 18 independent observations.** They are **9 matched pairs** — one high-risk and
one low-risk pool per (group × clutch) cell, processed together with the same reference gene. The
unit of analysis is the within-pair difference.

- Primary: **paired *t*-test across the 9 pairs**
- Robustness: **Wilcoxon signed-rank**
- Negative control: **sham pairs reported separately**
- Assumption check: Shapiro–Wilk on the within-pair differences

| Contrast | Pairs | Mean Δ (high − low) | 95% CI | *t* | *p* | Cohen's *dz* | Wilcoxon *p* |
|---|---|---|---|---|---|---|---|
| **All pairs (primary)** | 9 | +0.286 | [+0.008, +0.563] | *t*(8) = 2.375 | **0.045** | 0.79 | 0.055 |
| All pairs, log₂ scale | 9 | +0.350 | [+0.045, +0.656] | *t*(8) = 2.648 | **0.029** | 0.88 | 0.027 |
| Injured only | 6 | +0.391 | [−0.026, +0.807] | *t*(5) = 2.412 | 0.061 | 0.99 | 0.094 |
| **Sham only (control)** | 3 | +0.075 | [−0.334, +0.484] | *t*(2) = 0.790 | 0.512 | 0.46 | 0.500 |

Direction: 5 of 6 injured pairs run high > low.

**The sham control is what makes this interpretable.** If the pooling procedure, plate layout or
normalisation introduced a systematic high-vs-low difference, sham pairs would show it too. They do
not.

**Read the injured-only row honestly:** it has the largest effect size of the three but, at 6 pairs,
does not reach significance on its own. The all-pairs test is the primary one. The injured and sham
rows show *where* the effect sits — they are not two independent confirmations of it.

### 8.4 What was deliberately not done

**No regression of c-fos on a pooled continuous risk score across all 18 pools.** The risk score is
not on a comparable scale between dose groups — low and high dose move τ in opposite directions — so
pooling collapses the very contrast the pairing exists to isolate, and the relationship disappears.
Running it would produce a null result that means nothing about the biology.

### 8.5 Interpretation and limits

Larvae flagged high-risk by a purely behavioural model carry **27.5% more c-fos transcript** relative
to `rpl13a` than behaviourally low-risk siblings on the same plate (geometric mean ratio 1.275, 95%
CI [1.032, 1.575], back-transformed from the paired log₂ analysis — the appropriate scale for a fold
change). Elevated immediate early gene
expression in the absence of any provoking stimulus is what a chronically over-active network looks
like transcriptionally — the molecular counterpart of the reduced inhibitory gain that a long τ
reports behaviourally. **Two independent modalities, different fish, same latent variable.**

This is corroboration, not proof. It is a bulk measurement on pooled tissue: it cannot localise the
signal to a cell type or region, and it cannot distinguish loss of parvalbumin-positive interneuron
function from increased glutamatergic drive. Both are documented consequences of brain injury.

![Paired c-fos pools](results/figures/fig08_cfos_paired.png)

---

## 9. Secondary outcome: PTZ threshold — underpowered, and stated as such

Pentylenetetrazol is a non-competitive GABAₐ receptor antagonist that binds at the picrotoxin site,
reduces chloride conductance and removes inhibitory brake from the network. In larval zebrafish it
produces stereotyped, dose-dependent seizure behaviour with electrographic correlates (Baraban et
al., 2005), making it the standard pharmacological probe of seizure threshold. If injured larvae have
less inhibitory reserve, a fixed challenge dose should push more of them across it.

| Group | Seized / n | Proportion | Wilson 95% CI |
|---|---|---|---|
| Sham | 3/12 | 0.25 | [0.09, 0.53] |
| Low impact | 7/12 | 0.58 | [0.32, 0.81] |
| High impact | 7/10 | 0.70 | [0.40, 0.89] |

χ²(2) = 4.93, *p* = 0.085, Cramér's V = 0.381, n = 34.

> **This probe is underpowered and no conclusion in this project rests on it.** With 34 animals across
> three groups, post-hoc power for the observed sham-versus-injured difference (Cohen's *h* = 0.80) is
> **61%**. The direction is consistent with the primary result. That is all it is offered as.

---

## 10. Conclusion

The acute trajectory of startle-habituation kinetics — measured against each animal's own pre-injury
baseline and interpreted in the light of blast dose — separates injured larvae that develop
spontaneous burst activity from those that do not, at **AUC 0.833** and **76.5% total accuracy** under
leave-one-clutch-out nested cross-validation, with permutation *p* = 0.001. Sickness, injury severity,
and pre-existing trait were each tested as alternative explanations and each falls short. An
independent molecular assay on a separate cohort of animals points the same way.

What this supports is a mechanistic claim of modest scope: **a behavioural readout of inhibitory gain
in a defined sensorimotor circuit, sampled within 24 hours of injury, carries information about which
individuals are undergoing epileptogenesis.** It does not identify the cellular lesion, and it is one
experiment in one species at one developmental stage.

**Why that matters.** Anti-epileptogenic drug trials fail partly because they enrol everyone who was
injured, most of whom would never have developed epilepsy. A cheap, non-invasive, acute assay that
enriches a cohort for genuine epileptogenesis attacks that problem directly. The behavioural assay
costs **$0.09 and 0.9 operator-minutes per animal**.

## 11. Limitations

Stated plainly, because they bound the claim:

- **n = 81 with 36 events.** Events per variable = 9.0, at the accepted floor. Confidence intervals
  on the AUC and on every coefficient are wide and are reported rather than smoothed over.
- **Three clutches means three outer folds.** Fold-to-fold spread is estimated from three numbers;
  the SD across folds is indicative, not precise.
- **7 injured fish were dropped** for a missing post-injury session. This is a complete-case analysis;
  attrition (6.8% overall) is not obviously outcome-related but has not been modelled.
- **Conversion is a behavioural proxy.** Spontaneous burst activity at 6 dpf is not
  electrographically confirmed epilepsy. Field-potential recording from the optic tectum would be
  the confirmatory measurement.
- **The c-fos assay is bulk and pooled.** 9 pairs, 4 larvae per pool. It corroborates; it does not
  independently establish, and it cannot localise.
- **Cell-type resolution is absent.** The inhibitory-failure account in §4 is inference from circuit
  behaviour, not direct observation of interneuron loss in these animals.
- **PTZ is underpowered** and the study is not designed to support any claim from it.
- **One species, one age, one injury model.** Generalisation to mammalian TBI is a hypothesis.

## 12. Future work

1. **Electrophysiological confirmation** — tectal field recordings to establish that "conversion" is
   electrographic seizure activity.
2. **Cell-type resolution** — transgenic GABAergic reporter lines to test the interneuron-failure
   account directly, rather than inferring it from τ.
3. **Pharmacological rescue** — if τ indexes inhibitory gain, an agent that restores inhibition during
   the latent period should reduce conversion *and* renormalise τ. That is the experiment that would
   move this from biomarker to mechanism.
4. **External replication** — new clutches, ideally a second laboratory. Three clutches is the
   binding constraint on every uncertainty estimate here.
5. **Prospective enrichment test** — use the model to select a high-risk cohort *before* outcome is
   known, and check whether the observed conversion rate matches the predicted one.

## 13. ISEF compliance and ethics

- **Vertebrate animal research.** Under ISEF rules zebrafish are not considered vertebrate animals
  until 72 hours post-fertilisation. This study follows larvae to **6 dpf**, so it **is** regulated
  vertebrate animal research. It requires prior SRC/IACUC approval, Form 5B, a Qualified Scientist,
  and must be conducted at a Regulated Research Institution. **Verify current-year rules with your
  fair's SRC before starting — do not rely on this summary.**
- **Humane endpoints.** Blast injury and PTZ challenge both cause more than momentary distress and
  must be covered explicitly in the approved protocol, with anaesthesia and euthanasia methods
  specified.
- **Reduction.** The `cf_*` and `ptz_*` cohorts exist so that terminal assays do not consume animals
  from the followed cohort. Session logging (§ operational metrics) supports honest reporting of
  animal use: 133 larvae, 650 recording sessions, 9 animals lost.
- **Data integrity.** Every reported statistic is machine-generated into
  [`results/all_statistics.csv`](results/all_statistics.csv) by the pipeline. Nothing in RESULTS.md is
  hand-transcribed. Seed is fixed and reported; reruns are byte-identical.

---

## 14. Running it

```bash
pip install -r requirements.txt
python run_all.py
```

Options:

```bash
python run_all.py --permutations 200   # faster run (default 1000)
python run_all.py --skip-permutation   # skip the permutation test
python run_all.py --seed 123           # override the random seed
```

Full run ≈ 2 minutes, most of it the 1,000-iteration permutation test — each iteration reruns the
*entire* nested cross-validation.

### Data

The workbook is `realdata.xlsx` at the repository root, with seven sheets: `habituation_trials`,
`fish_features`, `outcomes_6dpf`, `cfos_cohort_features`, `cfos_pools`, `ptz_challenge`,
`session_log`. To point the pipeline elsewhere, change `DATA_XLSX` in
[`src/config.py`](src/config.py). `cfos_pools.delta_delta_ct` is aliased to `delta_ddct` on load, so
either column name works.

### Repository layout

```
run_all.py                  runner — executes every step in order
src/
  config.py                 seed, paths, design constants, plot style
  statsbook.py              append-only ledger; everything reported passes through it
  io_data.py                workbook loading and sheet validation
  habituation.py            STEP 1 — nonlinear curve fitting, agreement checks
  features.py               design matrix, model specifications
  modeling.py               ridge logistic, nested CV, permutation, bootstrap
  step2_prediction.py       STEP 2 driver — ablation table, leakage, coefficients
  step3_cfos.py             STEP 3 — paired c-fos pools (orthogonal validation)
  step4_descriptive.py      STEP 4 — descriptives, operations, PTZ
  plotting.py               all figures (300 dpi PNG)
  report.py                 assembles RESULTS.md from the ledger
results/
  figures/*.png             13 figures at 300 dpi
  tables/*.csv              intermediate tables (fits, folds, pairs, ...)
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
| `fig06` | Full-model coefficients with bootstrap CIs |
| `fig07` | Per-fold AUC by model |
| `fig08` | **9 paired c-fos lines, high vs low risk** |
| `fig09` | Habituation curves by group and timepoint, including t = −1 |
| `fig10` | τ by timepoint and group with error bars |
| `fig11` | Converter vs non-converter τ trajectories, split by dose |
| `fig12` | Operational metrics: throughput, operator time, cost, attrition |
| `fig13` | PTZ proportions with Wilson CIs |

### Reproducibility

Seed `20260809`, set in [`src/config.py`](src/config.py) and reported at the top of every run and in
RESULTS.md. Assumption checks print during the run tagged `[PASS]` or `[FLAG]` — including the two
that legitimately flag (non-normal baseline τ in the injured groups, backed by Kruskal–Wallis; one
of 650 fits at a parameter bound).

---

## Background literature

Selected reading behind the neurobiological claims above. These support the framing; they are not
citations of results generated here.

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

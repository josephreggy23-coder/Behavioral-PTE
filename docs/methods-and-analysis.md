# Methods and analysis record

This file records the analysis that the repository reproduces computationally. The author-supplied
acquisition protocol is now documented separately in
[experimental-methods.md](experimental-methods.md), together with its differences from the
canonical workbook. This is a post-project technical record, not a reconstructed pre-experiment
research plan.

## Reproducible computational pipeline

### 1. Habituation feature

For each fish-session, 30 trial distances are fit by bounded nonlinear least squares:

`distance(k) = A * exp(-(k - 1) / tau) + C`

Four starting values are attempted and the fit with the smallest sum of squared errors is retained.
The pipeline records convergence, R-squared, RMSE, parameter-bound status, and the standard error of
tau. It also runs a log-linear diagnostic to show why subtracting an estimated offset and logging the
remaining values is unsuitable for these data.

### 2. Prediction table

The primary set contains pressure-wave-exposed fish complete on the required predictors and outcome:

- `dose`: high impact = 1, low impact = 0;
- `pre_tau`: tau at -1 hour;
- `dtau_0.5`: tau at 0.5 hours minus pre-injury tau;
- `dtau_24`: tau at 24 hours minus pre-injury tau.

The delta features are centered and scaled within dose using statistics learned from the training
fold only, followed by standard scaling and L2-penalized logistic regression. This is
dose-conditioned preprocessing plus an additive dose term; it is not a formal interaction or
moderation model.

### 3. Validation

- Outer validation: three-fold `GroupKFold`, holding out one whole clutch at a time.
- Inner tuning: the two remaining clutches select the inverse penalty `C` from the fixed grid.
- Primary discrimination: pooled out-of-fold AUC, with mean and range of the three fold AUCs.
- Threshold metrics: accuracy, balanced accuracy, sensitivity, specificity, PPV, and NPV at 0.50.
- Conditional null: 1,000 label shuffles within clutch, rerunning nested validation each time.
- AUC interval: fish are resampled within each observed clutch while fixed out-of-fold predictions
  are reused. The interval is conditional on these clutches and omits model-refitting uncertainty.
- Coefficient intervals: fish are resampled within each observed clutch and the final fixed-penalty
  model is refit. These are not whole-cluster bootstrap intervals.

The pipeline compares the full model with restricted models. A restricted model cannot, by itself,
rule out every biological alternative represented by its variable; it only measures predictive
information available from that input under the specified validation.

### 4. Orthogonal molecular validation

`cfos_pools` is pivoted into nine supplied high/low pairs, one pair per group-clutch cell. Paired
tests are run on fold change and log2 fold change. This is the orthogonal validation layer because a
separate qPCR modality tests molecular concordance with the supplied risk strata. Its
precision is limited because pool labels and selection are not regenerated, raw Ct/QC data are
missing, three group-pairs share each clutch, and the injured-only contrast is imprecise.

### 5. PTZ comparison

The secondary PTZ table reports per-group proportions, Wilson intervals, a three-by-two Pearson
chi-square test, and a descriptive power calculation. With 34 fish, it is explicitly underpowered.

## Experimental documentation status

The newly supplied protocol adds AB strain, husbandry, developmental timing, infrared-rig hardware,
visual dark-flash timing, the general syringe pressure-wave procedure, tissue handling, qPCR
workflow, PTZ exposure timing, and stated randomization and blinding procedures. These additions
clarify that the author-supplied protocol identifies visual dark flashes rather than an acoustic
pulse.

The methods narrative is not a substitute for primary records, and several statements conflict with
the workbook. The analysis therefore follows the data-visible structure below.

| Component | Analysis-aligned interpretation | Still needed |
|---|---|---|
| Cohorts | 133 longitudinal, 86 molecular, and 34 PTZ IDs are non-overlapping | Linkage records if the cohorts were intended to overlap |
| Injury | Categorical low/high dose is available | Drop heights, mass, calibration, traces, pressure/impulse values |
| Habituation | Workbook has 30 trials at five timepoints; supplied protocol identifies visual dark flashes | Videos, tracking code, response window, `responded` threshold |
| Outcome | `converted` remains a supplied behavioral label | Exact velocity/burst algorithm and adjudication record; electrographic validation |
| qPCR | 18 four-larva pools form nine supplied high/low pairs | Raw Ct replicates, primers, efficiencies, QC and pool-selection rule |
| PTZ | Separate 34-animal, all-2.5-mM cohort | Vehicle records, scoring code, preparation, safety and disposition records |
| Oversight | All procedures are described as ending before 7 dpf | Research site, approval determination/identifier and welfare records |
| Contributions | AI use is recorded in `ai-use-disclosure.md` | Student, mentor, laboratory and software contribution statement |

The computational methods are reproducible from the repository. The experimental protocol is now
substantially better documented, but the experiment is not independently reproducible or
authenticated until the listed primary records are linked.

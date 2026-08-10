# Methods and analysis record

This file distinguishes what the repository can reproduce computationally from experimental methods
that still require primary documentation. It is a post-project technical record, not a reconstructed
pre-experiment research plan.

## Reproducible computational pipeline

### 1. Habituation feature

For each fish-session, 30 trial distances are fit by bounded nonlinear least squares:

`distance(k) = A * exp(-(k - 1) / tau) + C`

Four starting values are attempted and the fit with the smallest sum of squared errors is retained.
The pipeline records convergence, R-squared, RMSE, parameter-bound status, and the standard error of
tau. It also runs a log-linear diagnostic to show why subtracting an estimated offset and logging the
remaining values is unsuitable for these data.

### 2. Prediction table

The primary set contains blast-exposed fish complete on the required predictors and outcome:

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

### 4. Molecular comparison

`cfos_pools` is pivoted into nine supplied high/low pairs, one pair per group-clutch cell. Paired
tests are run on fold change and log2 fold change. This analysis is exploratory because pool labels
and selection are not regenerated, raw Ct/QC data are missing, three group-pairs share each clutch,
and the injured-only contrast is imprecise.

### 5. PTZ comparison

The secondary PTZ table reports per-group proportions, Wilson intervals, a three-by-two Pearson
chi-square test, and a descriptive power calculation. With 34 fish, it is explicitly underpowered.

## Experimental methods that are not documented

| Component | Required details |
|---|---|
| Animals | Strain, source, age at each procedure, sex if determinable, housing, husbandry, health monitoring, disposition |
| Design | Sample-size rationale, randomization, allocation concealment, blinding, inclusion/exclusion rules, experimental unit |
| Blast | Device, geometry, calibration, pressure/impulse units, sham exposure, anesthesia, recovery, humane endpoints |
| Startle assay | Stimulus waveform/intensity, inter-trial interval, camera/tracker, response window, `responded` threshold |
| Outcome | Recording method and duration, burst definition, blinded scoring, adjudication, electrographic status |
| qPCR | Tissue/collection time, extraction, primers, efficiencies, technical replicates, Ct exclusions, plate map, batch handling |
| PTZ | Preparation, exposure duration, observation/scoring protocol, censoring, safety, disposal, animal disposition |
| Contributions | What the student, mentor, laboratory, and software/AI each contributed |

Until those details are supplied, the computational methods are reproducible but the experiment is
not independently reproducible from this repository.

# Experimental methods and repository reconciliation

This record incorporates the experimental methods supplied by the project author on 11 August
2026. It documents the intended acquisition protocol and is separate from the executable analysis
record in [methods-and-analysis.md](methods-and-analysis.md). The repository does not contain the
raw videos, pressure traces, qPCR Ct files, allocation list, blinding key, dated protocol, or animal
oversight identifiers needed to authenticate these details independently.

Where the supplied narrative and canonical workbook disagree, this document states both versions.
All sample sizes, variables, and results reported by the pipeline follow the workbook and code.

## 1. Animals and husbandry

Wild-type AB zebrafish were maintained at 28.5 °C on a 14:10-hour light:dark cycle. Pairwise
crosses were set up the evening before collection, and embryos were collected within 30 minutes of
light onset to constrain developmental staging.

Fertilized embryos were raised at 50–60 per 100-mm dish in 30 mL E3 medium at 28.5 °C. Medium was
exchanged daily and nonviable embryos were removed. 1-phenyl-2-thiourea was not used because
recording was performed under infrared illumination and the protocol sought to avoid its reported
behavioral effects. At 3 days post-fertilization (dpf), larvae were screened for grossly normal
morphology and an intact touch-evoked escape response. They were then transferred individually to
24-well plates containing 1 mL E3 per well. Plate and well identifiers were intended to provide a
persistent animal identity. All procedures were reported to end before 7 dpf.

The supplied record does not identify the animal source or facility, routine health-monitoring
program, approval number, or final disposition records.

## 2. Experimental design and timing

The protocol targeted 135 larvae: 15 sham, 15 low-impact, and 15 high-impact animals in each of
three independent clutches. The canonical workbook contains 133 longitudinal animals at baseline
(45 sham, 45 low-impact, and 43 high-impact), so all longitudinal analyses use the 133 observed
animals. Clutches were collected on separate days. Clutch is the outer-validation group in the
prediction pipeline and the dependency unit in the molecular sensitivity analysis.

The supplied protocol places both a pre-injury session and pressure-wave exposure at 96 hours
post-fertilization, whereas the workbook codes baseline as −1 hour. The exact clock timing therefore
requires confirmation. The narrative lists post-injury sessions at 0.5 and 24 hours and a final
outcome recording at approximately 67 hours after injury (6.8 dpf). The workbook additionally
contains 1- and 5-hour sessions. It therefore represents five habituation timepoints: −1, 0.5, 1,
5, and 24 hours relative to injury. Across these timepoints, it contains 650 sessions and 19,500
individual trials.

The approximately 2.8-day interval between injury and outcome assessment was constrained by the
reported requirement to conclude procedures before 7 dpf.

## 3. Behavioral recording and tracking system

Recordings were acquired on a custom, top-down infrared video rig. The supplied protocol describes
the intended readout as whole-animal centroid displacement from each larva's silhouette. Because the
raw videos and tracking implementation are unavailable, the repository can verify and analyze the
derived trial distances but cannot reproduce the frame-to-centroid step.

- **Camera:** Raspberry Pi 4 Model B with a Raspberry Pi NoIR Camera Module v2 (Sony IMX219),
  mounted 240 mm above the plate. Acquisition was reported at 1280 × 720 pixels and 60 frames per
  second, providing approximately 10 pixels/mm across a 96-well-plate footprint.
- **Tracking illumination:** a 48-element, 850-nm infrared LED panel was positioned 80 mm below the
  plate. A 3-mm opal acrylic diffuser, 50 mm below the plate, homogenized the backlight.
- **Visual stimulus:** a 5000-K white LED panel was mounted 200 mm above the plate at 30° from
  vertical. Reported illuminance was 450 lux during lights-on periods and below 1 lux during each
  dark flash.
- **Switching and synchronization:** an IRLZ44N logic-level MOSFET connected to Raspberry Pi GPIO
  pin 18 switched the stimulus panel. A 3-mm red indicator LED in the camera field illuminated
  inversely with the panel to mark stimulus onset in the recorded frames.
- **Enclosure:** the rig was housed in a matte-black enclosure to suppress reflections and ambient
  light.

Python control software using `RPi.GPIO` and `picamera2` generated the stimulus sequence and a
timestamped session log. The supplied protocol reports timing jitter below 5 ms. Source acquisition
software, version identifiers, calibration images, and original session logs are not included in
the repository.

## 4. Pressure-wave injury

The supplied protocol describes the syringe pressure-wave method reported by Alyenbaawi et al. and
Kanyo et al. A single larva in 1 mL E3 was loaded into a 10-mL Luer-Lok syringe held vertically in a
rigid support. A calibrated weight was released through a guide tube onto the plunger. Low- and
high-impact conditions differed in drop height; weight mass, syringe volume, and medium volume were
held constant. Sham larvae underwent matched handling and residence time, with the released weight
stopped mechanically before impact.

The protocol also describes a piezoelectric pressure transducer mounted through the syringe wall
and sampled at 100 kHz, with peak overpressure and positive-pulse duration recorded for each event.
Those traces and per-animal pressure values are not present in the workbook. The current analysis
therefore uses the supplied categorical low-/high-impact group, encoded as a binary dose term; it
does not use measured peak overpressure.

Injury was reported to take approximately 45–60 seconds per larva. A priori exclusions were death,
visible cranial hemorrhage, absent touch response, or failure to right, applied while blinded to
later outcome. The repository does not contain event-level pressure logs, exclusion records, or the
pilot calibration data described in the narrative.

## 5. Visual dark-flash habituation assay

Larvae were placed in clear, flat-bottom 96-well plates containing 300 µL E3 per well. The protocol
specifies a 30-minute initial acclimation and 10 minutes before repeat sessions under 450-lux white
illumination. Each session comprised 10 minutes of baseline locomotor recording followed by 30
one-second periods of complete darkness separated by 15-second interstimulus intervals. This is a
**visual dark-flash assay**, not an acoustic startle assay.

For trial number `k`, response distance was modeled as

`response(k) = A × exp(−(k − 1) / tau) + C`

where `A` is response amplitude, `tau` is the habituation decay constant in trials, and `C` is the
asymptotic response floor. The supplied bounds were `A` in [0, 50], `tau` in [0.5, 40], and `C` in
[0, 20]. The executable pipeline instead uses wider bounds: `A` and `C` in [0, 200] and `tau` in
[0.1, 100]. It tries four starting values and retains the converged bounded nonlinear-least-squares
fit with the smallest residual sum of squares. Larger `tau` indicates slower habituation.

The workbook supplies trial distance and a binary `responded` field. The response-tracking window,
centroid-extraction implementation, and threshold used to create `responded` remain unavailable.

## 6. Later behavioral outcome

The protocol describes a one-hour, stimulus-free recording beginning approximately 67 hours after
injury, after 30 minutes of acclimation under dim illumination. High-velocity burst rate and mean
burst duration were retained as continuous outcomes. The narrative states that conversion was
prespecified as a burst rate above the 95th percentile of the sham distribution.

That rule does not reproduce the workbook's binary `converted` field: 9–11 of 133 labels disagree,
depending on the percentile convention. Consequently, the analysis continues to describe
`converted` as a supplied seizure-like behavioral label. The threshold implementation, velocity
cutoff, source videos, blinded scoring record, and any adjudication criteria are still needed before
the label can be regenerated.

## 7. Molecular cohort and qPCR

The supplied laboratory workflow states that plates were chilled immediately after outcome
recording to arrest transcription, followed by tricaine methanesulfonate overdose. Larvae were
transferred individually, snap-frozen in liquid nitrogen, and stored at −80 °C. Total RNA was
extracted with a silica-column kit and DNase I treatment, with an accepted A260/280 range of
1.9–2.1. Reverse transcription used a normalized 500-ng RNA input. SYBR Green qPCR targeted
`fosab`, used `rpl13a` as the reference gene, and included three technical replicates per pool.
No-template, no-reverse-transcriptase, melt-curve, amplification-efficiency, and reference-gene
stability controls were specified. Relative expression was defined by the 2^(−ΔΔCt) method.

The canonical workbook represents a non-overlapping molecular ID cohort of 86 `cf_*` animals, with
no linkage to the longitudinal IDs. Seventy-two of them appear in 18 pools of four larvae; 14 are
unpooled. This conflicts with the supplied narrative's estimate of three larvae per extraction and
approximately 27 pools. The repository therefore reports the workbook structure and treats the
nine supplied high-/low-risk pool pairs as the orthogonal molecular validation layer. The narrative
describes pooling by conversion status, while the workbook instead stores undocumented `high_risk`
and `low_risk` labels. Raw Ct replicates, primer sequences, efficiencies, QC outputs, and a
reproducible risk-pool selection rule are not available.

## 8. Pentylenetetrazol challenge

The supplied protocol specifies 2.5 mM pentylenetetrazol (PTZ) in E3, 30 minutes of locomotor
recording, latency to a high-velocity burst, and right-censoring at 1800 seconds. It describes three
washes in E3 after exposure and frames PTZ as a susceptibility probe rather than an
epileptogenesis model.

The workbook's 34 `ptz_*` identifiers do not overlap the longitudinal IDs and have no linkage to
them. The executable analysis therefore treats PTZ as a separate, secondary ID cohort—not as a
randomized arm within the 133-animal longitudinal cohort. It reports group proportions, Wilson
intervals, a three-by-two Pearson chi-square test, and a low-versus-high Fisher exact test. The
cohort is underpowered and no primary conclusion rests on it. Exposure preparation, vehicle
records, scoring code, safety documentation, and disposition records are not included.

## 9. Randomization, blinding, and reproducibility

The supplied narrative states that assignments were generated programmatically in advance,
interleaved across plate positions, and never allocated by column. Video files were reportedly
renamed to random codes before scoring, with the key unsealed only after scoring. It also states
that the conversion threshold and exclusions were documented before collection, and that operator
time, consumables, and attrition were logged by session.

The workbook includes session-level operational summaries, and the analysis holds out whole
clutches during validation. The allocation list, plate maps, blinded filenames, key, dated analysis
plan, and exclusion log are not present, so the randomization and blinding statements remain
author-supplied rather than independently auditable.

## 10. Computational analysis

The executable analysis is documented in [methods-and-analysis.md](methods-and-analysis.md) and
implemented in `src/`. It fits every 30-trial curve, constructs within-animal change features,
performs nested leave-one-clutch-out ridge-logistic validation, reruns that nested prediction
analysis under 1,000 within-clutch label permutations, and generates all tables, figures, and
reported statistics.

The supplied methods proposed additional analyses that are not implemented in the current code,
including a mixed-effects model of `tau`, continuous peak-overpressure prediction, a blocked linear
model for qPCR, and Kaplan–Meier/Cox modeling for PTZ. These should not be described as completed
unless the relevant source variables and code are added and the outputs regenerated.

## Reconciliation summary

| Item | Supplied narrative | Canonical workbook/current pipeline |
|---|---|---|
| Longitudinal sample | 135 targeted | 133 observed at baseline |
| Habituation schedule | −1, 0.5, and 24 h listed | −1, 0.5, 1, 5, and 24 h |
| Curve-fit bounds | A: 0–50; tau: 0.5–40; C: 0–20 | A/C: 0–200; tau: 0.1–100 |
| Injury predictor | Continuous peak overpressure | Categorical low/high dose only |
| Conversion rule | Burst rate above sham 95th percentile | Supplied label; 9–11/133 disagree depending on percentile convention |
| qPCR pooling | Three larvae/pool, approximately 27 pools, grouped by conversion | Four larvae/pool; 18 supplied high-/low-risk pools from 72 of 86 molecular IDs |
| PTZ design | Described as an assigned arm | Separate 34-animal cohort by identifier |
| Inferential methods | Mixed models, blocked qPCR model, survival analysis proposed | Ridge prediction, paired qPCR comparisons, PTZ contingency analyses implemented |

## References supplied with the protocol

1. Alyenbaawi H, et al. Seizures are a druggable mechanistic link between TBI and subsequent
   tauopathy. *eLife*. 2021;10:e58744.
2. Kanyo R, Wang CK, Locskai LF, Kim T, Allison WT. Delivering traumatic brain injury to larval
   zebrafish. *Methods in Molecular Biology*. 2023. PMID: 37668902.
3. A larval zebrafish model of traumatic brain injury: optimizing the dose of neurotrauma for
   discovery of treatments and aetiology. *Biology Open*. 2025;14(2):bio060601.
4. Beppi C, Penner M, Straumann D, Bögli SY. Biomechanical induction of mild brain trauma in larval
   zebrafish: effects on visual startle reflex habituation. *Brain Communications*.
   2023;5(2):fcad062.
5. Köcher L, Beppi C, Penner M, Meyer S, Bögli SY, Straumann D. Concussion leads to opposing
   sensorimotor effects of habituation deficit and fatigue in zebrafish larvae. *Brain
   Communications*. 2024;6(6):fcae407.
6. Marsden KC, Granato M. In vivo Ca²⁺ imaging reveals that decreased dendritic excitability drives
   startle habituation. *Cell Reports*. 2015;13:1733–1740.
7. Santistevan NJ, et al. Cacna2d3, a voltage-gated calcium channel subunit, functions in vertebrate
   habituation learning and the startle sensitivity threshold. *PLOS ONE*. 2022;17:e0270903.
8. Afrikanova T, et al. Validation of the zebrafish pentylenetetrazol seizure model: locomotor versus
   electrographic responses to antiepileptic drugs. *PLOS ONE*. 2013;8(1):e54166.
9. Cho SJ, et al. Zebrafish model of posttraumatic epilepsy. *Epilepsia*.
   2020;61(7):1197–1209.
10. Heylen L, et al. Pericardial injection of kainic acid induces a chronic epileptic state in
    larval zebrafish. *Frontiers in Molecular Neuroscience*. 2021.
11. Induction of seizures and initiation of epileptogenesis by pilocarpine in zebrafish larvae.
    *Frontiers in Molecular Neuroscience*. 2024.
12. Kumar A, Singh D. Temporal dynamics of c-Fos expression in a zebrafish larva model of
    pentylenetetrazole-induced seizures. 2025. doi:10.1177/15458547251401472.
13. Activity-dependent expression of neuronal PAS domain-containing protein 4 (`npas4a`) in the
    developing zebrafish brain. *Frontiers in Neuroanatomy*. 2014;8:148.
14. van Smeden M, et al. Sample size for binary logistic prediction models: beyond events per
    variable criteria. *Statistical Methods in Medical Research*. 2019.
15. Vabalas A, Gowen E, Poliakoff E, Casson AJ. Machine learning algorithm validation with a limited
    sample size. *PLOS ONE*. 2019;14(11):e0224365.

These references were supplied with the methods narrative and have not been independently checked
by this repository.

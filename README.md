# PSDM

Inference-only release of the final PSDM model described in the manuscript.
It contains one frozen ensemble of 350 LightGBM learners and the fitted
preprocessing parameters needed to apply it. Training, feature selection,
plotting and manuscript-analysis scripts are not included.

**Prepare annotations separately before prediction.** The VCF files contain
variant coordinates and reference labels, not model features. PSDM cannot
infer its sequence, conservation and splicing features from VCF coordinates
alone. No annotation matrices or stored prediction scores are bundled.

## Contents

```text
PSDM/
  README.md
  data/
    train.vcf          # 2,361 development variants
    test238.vcf        # 238 TestSet I variants
    external462.vcf    # 462 ClinVar evaluation variants
  code/
    PSDM.model.zip     # final model, preprocessing and feature schema
    predict.py         # model loading and prediction only
    requirements.txt
```

All VCF files use GRCh38, 1-based positions and biallelic SNVs. `LABEL=1`
and `LABEL=0` are the study's positive and negative reference classes;
labels are ignored by the predictor. The development VCF excludes the
single variant that failed the study's GRCh38 REF check.

## Run prediction

Use 64-bit Python 3.12. From the `PSDM` directory:

```bash
python -m pip install -r code/requirements.txt
python code/predict.py --describe > annotation_schema.json
python code/predict.py --vcf data/test238.vcf --annotations /path/to/my_annotations.tsv.gz --output output/test238_predictions.tsv
```

Replace `/path/to/my_annotations.tsv.gz` with **your own annotation file**.
The command does not download or generate annotations. Python dependencies
are required; no network access or GPU is needed for the prediction step.
The model is found relative to `predict.py`, so the script can also be
called from another working directory using absolute input paths.

Output columns are `variant_id`, `chrom`, `pos`, `ref`, `alt`, `input_id`,
`psdm_score` and `predicted_label`. The score is the mean of all 350 learner
probabilities; the positive rule is `score >= 0.5`. Existing output files
are not overwritten.

## Annotation input

The `--describe` command exports the complete feature schema and annotation
contract from the model. Supply a tab-separated table, optionally gzip
compressed, with the following two metadata lines:

```text
##reference=GRCh38
##PSDM_SCHEMA_SHA256=28ab2eea25779ccf536a0eaa742555d71d4104f58dd6248f5f0232cdf4c6d1b7
```

The next line must contain `chrom`, `pos`, `ref`, `alt` and the 260 feature
names listed by `--describe`. Supply **raw, unscaled values** under those
names; the fitted model performs imputation and scaling internally.
Extra columns are ignored. Column order and row order may differ: input
rows are matched to VCF records by chromosome, position, REF and ALT.
Both `1`/`chr1` style chromosome names are accepted. The VCF must declare
`##reference=GRCh38` or `##reference=hg38`; no genome-build conversion is made.

The 260 fields are the union used by the five folds: 41 iFeatureOmega,
44 ftrCOOL, 40 GPN-MSA and 135 fixed features. Each fold still uses its
original 231 inputs, so this reduction changes neither the fitted models
nor the features any learner receives. The fixed block contains 60 DNA
shape, 26 conservation and 49 RNA/splicing features.

Annotation definitions must match the study, including the sequence
windows and ref/alt/difference summaries, the GPN-MSA embedding operation
(not a scalar GPN variant score), CADD GRCh38-v1.6 conservation fields,
SpliceAI distance 50, and MMSplice 2.4.0 with GENCODE 43. Merely renaming
columns from another annotation pipeline does not establish compatibility.
The model's exported annotation contract provides these requirements.
The fixed-feature formulas and the `rna__0`–`rna__48` column meanings are
listed in the reference section below.

Use `NA` or `.` only for genuine unavailable annotation values. Missing
annotation rows, missing feature columns, duplicate variants, infinity,
multiallelic records, non-SNV alleles and entirely missing feature blocks
cause an error. Confirm the reference allele and synonymous consequence
in the upstream annotation step; the predictor does not bundle a reference
genome or consequence annotator. The schema checksum verifies the feature
contract, not the biological accuracy of a user-supplied annotation table.

## Load in Python

```python
import sys
sys.path.insert(0, "code")
from predict import PSDM, read_vcf

model = PSDM()
variants = read_vcf("data/test238.vcf")
features = model.read_annotations("my_annotations.tsv.gz", variants)
scores = model.predict(features, model.names)
```

## Reproduction scope

With the same raw annotations, this release reproduces final-model
predictions. Reference AUROCs are 0.947391 on TestSet I and 0.955102 on
ClinVar462. These are retrospective evaluations. Predicting `train.vcf`
with all 350 learners gives in-sample scores, **not** the paper's held-out
development OOF scores. This inference-only package does not reproduce
training, annotation generation or every manuscript figure and table.

The organization follows the separation between annotation and model
inference in the authors' [EECFS/CFDPSM](https://github.com/xialab-ahu/EECFS#5-use-cfdpsm-to-inference-other-ssnvs),
[SynMall/SynScore](https://github.com/xialab-ahu/SynMall/blob/4e9d31b4503423ac2167d45aaa2eec91287a1c43/Code/2-Whole-Genome-Prediction.py)
and [seDSM](https://github.com/xialab-ahu/seDSM/blob/005202745408eb1a2fe1de66c88916a16705055d/prediction.py)
repositories. Their implementations have not been copied into PSDM.



Code, model and data redistribution terms remain to be confirmed by the
authors before public release; this package does not assert a license or
an existing public repository.

## Fixed-feature definitions

These are input definitions, not stored annotation results. Keep the column
names exported by `--describe`; the semantic names below explain positional
RNA aliases. Derived values are calculated upstream, before the model's
fold-specific imputation or standardization.

### DNA shape and conservation

DNA shape names encode their shape, summary and allele (`ref`, `alt`,
`diff`). Summarize each DNAshapeR track over its available values in the
111-bp sequence. `sd` is the sample standard deviation (`ddof=1`);
`diff` is the alternate summary minus the corresponding reference summary.
Unavailable summaries remain missing.

The first ten conservation fields are the raw CADD GRCh38-v1.6 inclAnno
values named in the exported schema. For derived `phcons` summaries, use
`priPhCons,mamPhCons,verPhCons`; for `phylop`, use only
`priPhyloP,mamPhyloP,verPhyloP`; for `gerp`, use `GerpRS,GerpN,GerpS`.
For each group calculate mean, population standard deviation (`ddof=0`),
maximum, minimum and maximum-minus-minimum, ignoring missing values.
An all-missing group remains missing. `Conservation_mean_proxy` is the
missing-aware mean of `GerpRS,priPhyloP,mamPhyloP,verPhyloP`. Do not add
ZooPhyloP or siPhy fields to these calculations.

### RNA/splicing columns

`DS` in the derived formulas denotes the four SpliceAI delta scores
clipped below at zero, in AG, AL, DG, DL order. The first four raw input
columns remain unmodified. MMSplice modules are acceptorIntron, acceptor,
exon, donor and donorIntron; each module difference is alt minus ref.
All MMSplice raw fields must come from the same selected exon record.

Means and maxima skip missing values and remain missing if all inputs are
missing. Sums and comparison counts have the frozen builder's behavior:
an all-missing set contributes zero to those particular summaries.
Comparisons with missing values are false. Maximum-DS ties select the
first entry in AG, AL, DG, DL order; if all DS values are missing, the
position selector uses AG. `log1p(x)` means `ln(1+x)`. These rules describe
genuine annotation missingness; a failed annotation job is an error.

| Input column | Meaning | Value or formula |
|---|---|---|
| `rna__0` | `DS_AG` | Raw DS_AG |
| `rna__1` | `DS_AL` | Raw DS_AL |
| `rna__2` | `DS_DG` | Raw DS_DG |
| `rna__3` | `DS_DL` | Raw DS_DL |
| `rna__4` | `DP_AG` | Raw DP_AG |
| `rna__5` | `DP_AL` | Raw DP_AL |
| `rna__6` | `DP_DG` | Raw DP_DG |
| `rna__7` | `DP_DL` | Raw DP_DL |
| `rna__8` | `MMS_delta_logit_psi` | Raw MMS_delta_logit_psi |
| `rna__9` | `MMS_ref_acceptorIntron` | Raw MMS_ref_acceptorIntron |
| `rna__10` | `MMS_alt_acceptorIntron` | Raw MMS_alt_acceptorIntron |
| `rna__11` | `MMS_ref_acceptor` | Raw MMS_ref_acceptor |
| `rna__12` | `MMS_alt_acceptor` | Raw MMS_alt_acceptor |
| `rna__13` | `MMS_ref_exon` | Raw MMS_ref_exon |
| `rna__14` | `MMS_alt_exon` | Raw MMS_alt_exon |
| `rna__15` | `MMS_ref_donor` | Raw MMS_ref_donor |
| `rna__16` | `MMS_alt_donor` | Raw MMS_alt_donor |
| `rna__17` | `MMS_ref_donorIntron` | Raw MMS_ref_donorIntron |
| `rna__18` | `MMS_alt_donorIntron` | Raw MMS_alt_donorIntron |
| `rna__19` | `MMS_pathogenicity` | Raw MMS_pathogenicity |
| `rna__20` | `MMS_efficiency` | Raw MMS_efficiency |
| `rna__21` | `MMS_diff_acceptorIntron` | MMS_alt_acceptorIntron - MMS_ref_acceptorIntron |
| `rna__22` | `MMS_diff_acceptor` | MMS_alt_acceptor - MMS_ref_acceptor |
| `rna__23` | `MMS_diff_exon` | MMS_alt_exon - MMS_ref_exon |
| `rna__24` | `MMS_diff_donor` | MMS_alt_donor - MMS_ref_donor |
| `rna__25` | `MMS_diff_donorIntron` | MMS_alt_donorIntron - MMS_ref_donorIntron |
| `rna__26` | `SpliceAI_DS_max` | max of clipped DS_AG,DS_AL,DS_DG,DS_DL, ignoring NaN |
| `rna__27` | `SpliceAI_DS_sum` | sum of clipped DS values, skip NaN; all missing -> 0 |
| `rna__28` | `SpliceAI_gain_max` | max of clipped DS_AG and DS_DG, skip NaN |
| `rna__29` | `SpliceAI_loss_max` | max of clipped DS_AL and DS_DL, skip NaN |
| `rna__30` | `SpliceAI_n_nonzero_001` | count of clipped DS > 0.01; NaN comparisons false |
| `rna__31` | `SpliceAI_n_high_02` | count of clipped DS >= 0.2; NaN comparisons false |
| `rna__32` | `SpliceAI_n_high_05` | count of clipped DS >= 0.5; NaN comparisons false |
| `rna__33` | `SpliceAI_abs_DP_at_max_DS` | abs(DP) paired with largest clipped DS in order AG,AL,DG,DL; ties choose first; all DS missing selects AG after NaN -> -infinity |
| `rna__34` | `SpliceAI_near_boundary_at_max_DS` | float(abs_DP_at_max_DS <= 10); NaN comparison false -> 0 |
| `rna__35` | `SpliceAI_log_abs_DP_at_max_DS` | log1p(abs_DP_at_max_DS) |
| `rna__36` | `SpliceAI_DSmax_x_near_boundary` | SpliceAI_DS_max * SpliceAI_near_boundary_at_max_DS |
| `rna__37` | `SpliceAI_DSmax_over_logDP` | SpliceAI_DS_max / (1 + SpliceAI_log_abs_DP_at_max_DS) |
| `rna__38` | `MMSplice_ref_mean` | mean of 5 MMS_ref module scores, skip NaN |
| `rna__39` | `MMSplice_alt_mean` | mean of 5 MMS_alt module scores, skip NaN |
| `rna__40` | `MMSplice_diff_mean` | mean of 5 MMS_diff module scores, skip NaN |
| `rna__41` | `MMSplice_abs_diff_max` | max(abs(5 MMS_diff module scores)), skip NaN |
| `rna__42` | `MMSplice_abs_diff_sum` | sum(abs(5 MMS_diff module scores)), skip NaN; all missing -> 0 |
| `rna__43` | `MMSplice_alt_minus_ref_mean` | MMSplice_alt_mean - MMSplice_ref_mean |
| `rna__44` | `MMSplice_abs_delta_logit_psi` | abs(MMS_delta_logit_psi) |
| `rna__45` | `MMSplice_abs_efficiency` | abs(MMS_efficiency) |
| `rna__46` | `MMSplice_path_x_abs_delta` | MMS_pathogenicity * MMSplice_abs_delta_logit_psi |
| `rna__47` | `SpliceAI_DSmax_x_MMS_pathogenicity` | SpliceAI_DS_max * MMS_pathogenicity |
| `rna__48` | `SpliceAI_DSsum_x_MMS_absdiff` | SpliceAI_DS_sum * MMSplice_abs_diff_sum |

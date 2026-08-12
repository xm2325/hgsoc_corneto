# GSE208216 organoid and DepMap validation contract

## Purpose

This work package asks whether Taylor OCM candidate genes or network features transfer to
independent HGSOC culture models. It does **not** treat all ovarian cell lines as HGSOC and it does
not treat DepMap dependency as proof of patient-level therapeutic efficacy.

## GSE208216 evidence boundary

The official GEO record exposes one raw-count matrix with 14 columns: 11 HGSOC patient-derived
organoids (PDO1, PDO2, PDO3, PDO5, PDO6, PDO7, PDO8, PDO10, PDO11, PDO12, and PDO15) and three
fallopian-tube organoids (FT1, FT8, and FT9). The archive checksum, exact sample set, integer-count
schema, and 33,610 gene rows are fixed in `config/external_validation/gse208216.json`.

The article reports 15 continuous PDO models, but that number describes the wider model resource,
not the deposited bulk-RNA matrix. Numeric copy-number-signature activities are not in the GEO
expression deposit. The paper points to controlled-access EGA study EGAS00001007189 for sWGS and
scDNA data. Consequently:

- the public expression analysis may compare the 11 deposited PDOs with the three FT organoids;
- paper-explicit qualitative labels can be used with an exact citation;
- a complete numeric `sample × copy-number-signature` table must remain unavailable until the
  controlled data or an author-supplied source table is obtained;
- values must not be digitized from a stacked-bar figure and presented as source measurements.

## DepMap positive set and comparators

`config/external_validation/depmap_hgsoc_models.tsv` defines a deliberately narrow positive set:
KURAMOCHI, OVSAHO, COV362, OVCAR4, and SNU119. Every run must resolve the stored ACH identifier to
the expected model name in that release's `Model.csv`. Missing or changed mappings are fatal.

Other ovary/fallopian-tube models form an explicitly labelled comparator group,
`other_ovarian_not_hgsoc_positive`. They are not silently promoted to HGSOC on the basis of lineage,
an ovarian collection site, or a generic serous label. In particular, commonly used ovarian models
such as SKOV3 or A2780 must not be used as positive HGSOC baselines without independent evidence.

DepMap releases are mutable. Record the release label, source URLs, file sizes, and SHA-256 hashes
for `Model.csv` and `CRISPRGeneEffect.csv` in every receipt. Do not combine files from different
releases.

The executable preflight intentionally stores only DepMap's official download landing page. If an
explicit quarterly release, `Model.csv`, `CRISPRGeneEffect.csv`, and that release's README are not
all present, it emits `status=blocked` with `scientific_success=false`. It never guesses a portal
asset URL or silently substitutes an older release.

On Roihu, `hpc/roihu/external_gse208216_fetch_audit.sbatch` atomically downloads and checksum-gates
the public GEO matrix. `hpc/roihu/external_depmap_preflight.sbatch` writes the explicit ready/blocked
DepMap receipt from operator-supplied environment variables. Neither script submits successor jobs.

## Analysis contract

1. Freeze a Taylor-derived `gene_symbol` candidate table before viewing external effects.
2. Audit the independent input file and emit hashes and model/sample identities.
3. For GSE208216, use an explicit, hashed Ensembl-to-symbol table; produce model-level
   `log2(CPM + 1)` values and descriptive PDO-versus-FT summaries. Do not infer symbols from row
   order or query a mutable service during the scientific run.
4. Extract candidate gene effects for the curated positive set and separately for other ovarian
   models. Missing candidates and missing screened models remain explicit.
5. Report per-gene model-level values and descriptive group summaries. With only five curated
   models, prioritize effect sizes and leave-one-model-out stability over thresholded p-values.
6. Interpret a more negative Chronos gene effect as stronger in-vitro dependency only. Concordance
   supports model-level functional plausibility; discordance falsifies broad transfer but can also
   reflect expression/model/platform differences.
7. Do not claim drug sensitivity, clinical response, or causality from gene effect alone.

The public DepMap download page states that `Model.csv` maps release metadata by `ModelID`, while
`CRISPRGeneEffect.csv` contains post-Chronos model-by-gene effects. File names and schemas must be
verified against the downloaded release README rather than assumed from an older quarter.

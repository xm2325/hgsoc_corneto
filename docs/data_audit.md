# Public-data audit

Snapshot date: 2026-08-11.

Generated TSV files use the literal `NA` for missing scalar values. Empty cells
in byte-exact upstream MAGE-TAB snapshots are preserved as received.

## RNA archives

| Study | ENA project | Runs | FASTQs | Compressed bytes | GiB |
|---|---|---:|---:|---:|---:|
| E-MTAB-7223 | PRJEB28709 | 36 | 72 | 188,841,300,121 | 175.872 |
| E-MTAB-10801 | PRJEB46736 | 36 | 72 | 172,318,782,780 | 160.484 |
| E-MTAB-11000 | PRJEB47842 | 12 | 24 | 33,440,119,254 | 31.144 |
| E-MTAB-14568 | PRJEB81794 | 33 | 66 | 83,162,442,959 | 77.451 |
| **Total** | | **117** | **234** | **477,762,645,114** | **444.951** |

Counts and byte totals are recalculated from the committed ENA `read_run`
reports, not copied from prose.

## Tighe phenotype boundary

Public supplementary assets contain:

- Table S1: 83 OCM rows from 68 patients, including histotype, biopsy type,
  chemotherapy-naive status, and TP53 annotations.
- Table S2: normalized ABCB1 read count for the same OCM panel.

They do **not** contain a machine-readable table of exact OCM-level paclitaxel
AUC/GI50, carboplatin AUC/GI50, or cumulative paclitaxel exposure. Values in
figures are not digitized into the primary analysis.

## Meeson input boundary

The pinned `katemeeson/PhD_2024` repository contains:

- the public sequential integration implementation;
- a notebook that restores an untracked `ocm_t` variable and then reads a local
  `bc_rnaseq_with_growths_adjusted.csv` path;
- a 49-row `OCM_clusters_and_fluxes.csv` output table; and
- cell-line dependency/TPI1 validation outputs.

The exact 49-OCM expression-plus-growth input is not committed upstream.
Consequently, exact numerical regeneration of those 49 models is distinguished
from (a) auditing the published sequential algorithm on fully public inputs and
(b) rebuilding OCM inputs from the four archived RNA studies.

## Meeson Human-GEM reconstruction

The Meeson paper reports that `Human-GEM-annotated.xml` was accessed in
September 2020 and contains 13,096 reactions and 3,628 genes. That locally named
file is not committed in any of the three public Meeson repositories. The
underlying public model can nevertheless be identified without guessing:

- Human-GEM `v1.4.1` was released on 2020-07-29 and was the latest release in
  September 2020.
- Its SBML has SHA-256
  `57d1b137f0c90d83a3e4f9a8225d74d37523594e6ee99f622b160a014d9f7050`.
- It has exactly 13,096 reactions and 3,628 genes.
- Its GPR partition is exactly 653 AND, 129 mixed AND/OR, 3,972 OR, 3,282
  single-gene, and 5,060 no-gene reactions, matching the public notebook.
- COBRApy decodes the SBML reaction identifier `R_biomass_human` to
  `biomass_human`, also matching the notebook.

Human-GEM `v1.4.0` has the same biological fingerprint and its XML differs from
`v1.4.1` only in metadata timestamps. We pin `v1.4.1` because it is the
chronologically correct public release for the reported access date.

The public integration code remains order-dependent by construction: it tests
candidate reaction bounds one at a time, retains a candidate only when the
current model remains above the growth threshold, and otherwise reopens that
reaction before moving to the next. Mixed AND/OR rules are classified but never
integrated. Reopened bounds are hard-coded to 0/±1000 rather than restored from
the input model. Our benchmark therefore reports both the published semantics
and a bounds-safe audit semantics, and never silently treats them as equivalent.

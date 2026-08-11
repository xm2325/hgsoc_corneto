# Regulatory CORNETO pilot

This track is a response-blind, reproducible network-inference pilot. It uses
public signed interaction priors, not a new phenotype table and not a causal
Taxol claim.

## Public inputs

`scripts/fetch_regulatory_sources.py` downloads the human CollecTRI and
OmniPath tables from their public OmniPath endpoints. It filters to directed,
unambiguous consensus stimulation/inhibition edges, removes self-loops and
preserves downloaded TSVs and normalizes source/target/sign at read time. The receipt records the
retrieval timestamp, URL, response metadata, SHA256, row counts and skipped-row
counts. Raw and normalized tables belong on Roihu project scratch; they are
not committed to GitHub. Source-specific licensing and attribution terms must
be retained with any redistribution.

## Inference policy

The pilot computes signed CollecTRI regulon z-scores from log1p TPM, selects
deterministic extreme TF outputs, chooses bounded upstream expression-derived
priors, and emits normalized signed-PKN edges; the CORNETO/CARNIVAL solve is explicitly blocked.
Expression-derived priors are not perturbation experiments. Each sample records
the source hashes, selected nodes/edges, corneto blocked status and any blocked/error
reason. A blocked sample is not silently treated as a network result.

The first pilot is E-MTAB-14568 primary tumour samples (8–12 samples). Only
after the source receipt and small pilot are stable should the same frozen
policy be expanded to the 60-sample primary tumour cohort. Cross-sectional,
longitudinal and acquired-resistance interpretations remain descriptive until
the exact phenotype tables and patient-grouped validation are available.

Official method references: [CORNETO CARNIVAL](https://corneto.org/latest/guide/signaling/carnival.html),
[CollecTRI](https://github.com/saezlab/collecTRI), and the
[OmniPath API](https://omnipath.readthedocs.io/en/latest/api.html).

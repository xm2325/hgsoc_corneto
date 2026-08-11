# Methodological boundaries

## Track A

Drug response is downstream evaluation data. It is not available to CollecTRI,
PKN construction, CORNETO optimization, hyperparameter stability analysis, or
module definition. Any supervised tuning is nested inside patient-grouped
training folds.

## Track B

The published Meeson algorithm sequentially proposes expression-derived
reaction constraints and retains a constraint only when the model remains above
an experimental growth threshold. Reaction order can therefore matter in
principle and must be tested empirically.

The proposed binary retain/relax formulation optimizes constraint retention
globally, subject to steady state and sample-specific growth feasibility. The
joint extension adds a structured penalty over reaction-selection indicators
across samples. This is new work. It must not be described as the exact CORNETO
MitoCore sFBA implementation or as an exact reproduction of Meeson unless the
same input matrix, model revision, media, solver, and numerical tolerances are
available.

## Reporting

Primary stability objects are recurrent modules, reaction/gene selection
frequencies, and solution ensembles. One arbitrary optimal flux vector or one
selected signalling edge is insufficient evidence of a stable mechanism.

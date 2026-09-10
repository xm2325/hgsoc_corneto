# Revised model validation, 2026-09-10

User authorization: revise the model and revalidate; no old b25 hold is released.

This first revision tests a **GPR-constrained parsimonious FBA** alternative,
not minimum-cardinality CORNETO. It deliberately removes binary-indicator
acceptance from this diagnostic: a reaction is selected from the actual saved
flux (absolute threshold 1e-7), then every unselected reaction is strictly zeroed
and an independent continuous LP must recover the specified growth floor.
This avoids claiming that a numerically inconsistent indicator vector is a network.
It does not by itself fix or validate the old binary formulation.

Frozen diagnostic panel: the same four response-blind OCMs in the audited r3
smoke receipt, one per cohort. The runner records source/model/helper hashes,
selected reactions, fluxes, residuals and the independent support LP result.
GPR AND=min, OR=sum; positive log1p(TPM) caps at scale1 are retained as an
explicit uncalibrated assumption. For each sample, first maximize biomass;
at 50% and 90% of that sample's maximum, minimize total absolute flux.
These are prespecified growth sensitivities, not fitted to experimental growth.

All eight support cases must pass: optimal LPs, mass/bound residual <=1e-6,
fixed-support maximum >= intended growth floor minus1e-6. Four negative controls
close uptake in every one-metabolite boundary; positive growth there blocks
completion and requires boundary/stoichiometric diagnosis, not automatic repair.
Each LP has a30-second limit; the job uses4CPUs/32G/15min, no Gurobi session.
Outputs are fresh and every case is persisted. No canonical downstream release.

## Medium limitation and next gate

The Nelson/Taylor paper establishes OCMI use, but the current repository's recipe
anchor is explicitly incomplete. [Primary source](https://www.nature.com/articles/s41467-020-14551-2).
This pilot retains legacy bounds for a numerical comparison, not a reconstructed
OCMI medium. Closing three energy uptakes previously failed to change growth;
that result must not be used to pronounce the medium valid. Before biological
or full-cohort release, complete basal-product/chemical mapping and explicit
serum assumptions, then assess uptake-rate sensitivity. Concentrations are not
fluxes. The user-requested model revision therefore remains staged until this
medium gate and any subsequent binary CORNETO formulation tests pass.

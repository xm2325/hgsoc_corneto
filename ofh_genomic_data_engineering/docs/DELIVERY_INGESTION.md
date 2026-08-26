# Provider delivery ingestion contract

## Purpose

The genomic pipeline should not treat a downloaded object as valid merely because it is readable. A provider delivery is admitted only after its manifest and delivered VCF agree on the fields the pipeline can verify directly.

## Manifest

`config/delivery_manifest.json` contains:

- `delivery_id`: stable identifier for one provider delivery;
- `provider`: source organisation or delivery producer;
- `source_uri`: URI expected by the pipeline;
- `source_format`: currently `vcf.gz`;
- `reference_genome`: provider-declared assembly;
- `sha256`: expected bytes-level SHA-256 of the delivered object;
- `sample_roster.count`: expected number of VCF samples;
- `sample_roster.ids_sha256`: SHA-256 of sorted sample IDs joined with newlines.

The roster hash makes sample-set checking compact without storing participant identifiers again in a second file.

## Admission decisions

`validate_delivery.py` writes a structured JSON decision.

- `PROCESS`: manifest/source checks passed and the delivery ID is not present in the supplied registry.
- `NOOP`: an optional registry already contains the same delivery ID, source SHA-256 and delivery fingerprint. `should_process` is false.
- `REJECT`: any required contract fails, or the same delivery ID is registered with different content.

The default public CI run has no persistent registry, so its real-data delivery is expected to return `PROCESS`. Registry behaviour is covered with deterministic unit tests.

## Fail-closed checks

The validator rejects:

- missing or malformed required manifest fields;
- unexpected source URI when the pipeline supplies one;
- unsupported source format;
- declared genome build that does not match pipeline policy;
- empty source files;
- source SHA-256 mismatch;
- unreadable or malformed VCF headers;
- duplicate VCF sample IDs;
- sample-count mismatch;
- sample-roster hash mismatch;
- duplicate registry entries for one delivery ID;
- delivery-ID reuse with changed content.

## Release binding

A passing delivery receives a deterministic `delivery_fingerprint` computed from the validated delivery identity, source hash, reference-genome declaration and observed sample-roster evidence.

The provenance stage refuses failed/non-`PROCESS` decisions, checks the validated source SHA against the downloaded-source SHA, then records the delivery fingerprint. The fingerprint is also included in provenance parameters, which are part of the deterministic release-ID basis. The release stage checks the delivery fingerprint and source SHA against provenance before the normal release validator runs.

## Boundary

The `reference_genome` field is a provider declaration checked against pipeline policy. This validator does not infer GRCh37 versus GRCh38 cryptographically from variant coordinates.

The optional registry file is an interface for idempotency semantics, not a persistent multi-user registry service. A production deployment would normally back the same decision contract with durable metadata storage and concurrency control.

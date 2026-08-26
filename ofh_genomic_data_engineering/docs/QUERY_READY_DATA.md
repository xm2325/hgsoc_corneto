# Query-ready genomic data layer

## Purpose

The pipeline does not treat Parquet as a file-extension conversion. The v0.3 data layer assigns stable physical types, records the Parquet schema, uses ZSTD compression and validates region filtering with a second query engine.

## Type contract

Identifiers and chromosome labels remain strings. Genomic positions and count fields are stored as `int64`. Frequencies, missingness fractions, Hardy-Weinberg probabilities and principal-component scores are stored as floating-point values.

The type conversion is strict: a non-numeric value in a required numeric PLINK field stops the export instead of silently becoming a string. The generated `schema_manifest.json` records each table's row count and Arrow data types.

The variants table is sorted by chromosome, genomic position, reference allele and alternate allele before Parquet output. This provides stable ordering and supports Parquet statistics for range filtering. No latency or scale claim is made from this small public fixture.

## Independent query contract

`query_contract.py` reads `variants.parquet` with DuckDB and checks a deterministic one-megabase range from the first observed genomic position. It verifies that:

- `POS` is physically stored as `int64`;
- the dataset contains variants and a valid genomic position range;
- DuckDB and pandas return the same count for the region predicate;
- returned genomic positions are ordered when requested.

The process writes `query_validation.json` and exits non-zero on failure. The release process depends on this successful query result, so a candidate release is not accepted when the analysis data layer cannot satisfy its query contract.

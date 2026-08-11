# Roihu execution policy

All dependency-heavy metabolic and RNA work runs on the Roihu CPU partition.
The local checkout is used for source editing and lightweight metadata tests
only.

The frozen deployment location is:

```text
/scratch/project_2012997/xiaomei/hgsoc_corneto
```

The pinned public Human-GEM and Meeson archives are stored separately at:

```text
/scratch/project_2012997/xiaomei/hgsoc_corneto_inputs
```

This prevents upstream data from being accidentally committed or copied into
multiple analysis directories.

## Environment and smoke benchmark

From `roihu-cpu.csc.fi`:

```bash
bash hpc/roihu/setup_environment.sh
sbatch hpc/roihu/meeson_smoke.sbatch
```

The setup uses CSC's `python-data/3.12` module and a CPU-specific virtual
environment. CORNETO is installed from the immutable commit in
`config/sources.yaml`; Human-GEM and the Meeson repository archive are verified
by SHA-256 before use.

Job monitoring is deliberately stage-based. Record the job ID at submission,
then check once near the requested wall-time or after Slurm reports a terminal
state. Do not poll `squeue` or `sacct` in a tight loop. Inspect the single job log
and the generated completion receipt together.

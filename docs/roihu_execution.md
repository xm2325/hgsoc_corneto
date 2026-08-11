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

## Environment and solver smoke benchmark

From `roihu-cpu.csc.fi`:

```bash
bash hpc/roihu/setup_corneto_environment.sh
```

The setup uses CSC's `python-data/3.12`, the dedicated
`/scratch/project_2012997/xiaomei/hgsoc_corneto_env` environment, the pinned
CORNETO commit in `config/sources.yaml`, and Human-GEM/Meeson archives
verified by SHA-256 before use. Commercial solver clients are optional:

```bash
export GRB_LICENSE_FILE=/scratch/project_2012997/xiaomei/hgsoc_corneto/gurobi_xm.lic
INSTALL_GUROBI=1 bash hpc/roihu/setup_corneto_environment.sh
```

The license file remains on Roihu and is ignored by Git; receipts record only
package/license status, never license content. A tiny compute-node smoke should
pass before a real pilot:

```bash
sbatch --export=ALL,SOLVER=gurobi,GRB_LICENSE_FILE \
  hpc/roihu/gurobi_license_smoke.sbatch
```

The 14568 pilot is response-blind metabolic CORNETO plumbing, not a final
Taxol/regulatory result:

```bash
sbatch --export=ALL,SOLVER=gurobi,GRB_LICENSE_FILE \
  hpc/roihu/corneto_14568_pilot.sbatch
```

Job monitoring is deliberately stage-based. Record the job ID at submission,
perform a startup gate after 30--60 seconds, then check at planned milestones
or after Slurm reports a terminal state. Do not poll `squeue` or `sacct` in a
tight loop. Inspect the single job log and the generated completion receipt
together.

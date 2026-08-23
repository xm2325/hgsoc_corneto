#!/usr/bin/env python3
"""Run time-bounded checkpoint solves and always preserve solver telemetry."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
from typing import Any


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected a JSON object: {path}")
    return value


def _atomic_write(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def _context(path: Path) -> tuple[dict[str, Any], str]:
    value = _read_json(path)
    if value.get("status") != "prepared":
        raise ValueError("scientific context is not prepared")
    return value, _sha256(path)


def _bounds(
    context: dict[str, Any], conditions: list[str]
) -> dict[str, dict[str, tuple[float, float]]]:
    return {
        condition: {
            reaction: (float(interval[0]), float(interval[1]))
            for reaction, interval in context["reaction_bounds"][condition].items()
        }
        for condition in conditions
    }


def _completed_matches(path: Path, context_sha: str, condition: str | None = None) -> bool:
    if not path.is_file():
        return False
    value = _read_json(path)
    expected = {"status": "completed", "context_sha256": context_sha}
    if condition is not None:
        expected["condition"] = condition
    if all(value.get(key) == wanted for key, wanted in expected.items()):
        return True
    raise ValueError(f"existing canonical receipt is incompatible: {path}")


def _attempt_path(directory: Path, stem: str) -> Path:
    return directory / f"{stem}{_attempt_suffix()}.json"


def _attempt_suffix() -> str:
    job = os.environ.get("SLURM_JOB_ID", "local")
    task = os.environ.get("SLURM_ARRAY_TASK_ID")
    return f"_job{job}" + (f"_task{task}" if task is not None else "")


def independent(args: argparse.Namespace) -> int:
    context, context_sha = _context(args.context)
    conditions = list(context["conditions"])
    if args.array_index < 0 or args.array_index >= len(conditions):
        raise IndexError(f"array index {args.array_index} outside 0..{len(conditions) - 1}")
    condition = conditions[args.array_index]
    canonical = args.output_dir / f"{args.array_index:03d}_{condition}.json"
    if _completed_matches(canonical, context_sha, condition):
        print(json.dumps({"status": "existing_valid", "output": str(canonical)}))
        return 0

    import cobra

    from hgsoc_corneto.metabolic.instrumented_fba import solve_independent_instrumented

    model = cobra.io.read_sbml_model(str(args.human_gem))
    if _sha256(args.human_gem) != context["model"]["sha256"]:
        raise ValueError("Human-GEM SHA256 differs from the frozen context")
    stem = f"{args.array_index:03d}_{condition}"
    artifact_prefix = args.artifact_dir / f"{stem}{_attempt_suffix()}"
    result = solve_independent_instrumented(
        model,
        condition=condition,
        objective=context["objectives"][condition],
        reaction_bounds=_bounds(context, [condition])[condition],
        independent_lambda=float(context["objective"]["independent_lambda"]),
        max_seconds=args.max_seconds,
        mip_gap=args.mip_gap,
        threads=args.threads,
        artifact_prefix=artifact_prefix,
        log_file=artifact_prefix.with_suffix(".gurobi.log"),
    )
    payload = {
        "status": result.status,
        "scientific_success": result.scientific_success,
        "schema_version": "metabolic_independent_checkpoint.v1",
        "condition": condition,
        "context_sha256": context_sha,
        "study_accession": context["study_accession"],
        "array_index": args.array_index,
        "solver": "gurobi",
        "independent_lambda": context["objective"]["independent_lambda"],
        "solution": (result.to_dict()["summaries"][0] if result.summaries else None),
        "instrumentation": result.to_dict(),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
        "claim_limit": (
            "Canonical result only when status=completed; partial incumbent is optimization "
            "diagnostic evidence and not a biological result."
        ),
    }
    if result.scientific_success:
        _atomic_write(canonical, payload)
        output = canonical
        exit_code = 0
    else:
        output = _attempt_path(args.attempt_dir, stem)
        _atomic_write(output, payload)
        exit_code = 2
    print(json.dumps({"status": result.status, "output": str(output)}))
    return exit_code


def joint(args: argparse.Namespace) -> int:
    context, context_sha = _context(args.context)
    if _completed_matches(args.output, context_sha):
        print(json.dumps({"status": "existing_valid", "output": str(args.output)}))
        return 0
    conditions = list(context["conditions"])

    import cobra

    from hgsoc_corneto.metabolic.instrumented_fba import solve_joint_instrumented

    model = cobra.io.read_sbml_model(str(args.human_gem))
    if _sha256(args.human_gem) != context["model"]["sha256"]:
        raise ValueError("Human-GEM SHA256 differs from the frozen context")
    artifact_prefix = args.artifact_dir / f"joint{_attempt_suffix()}"
    result = solve_joint_instrumented(
        model,
        objectives={condition: context["objectives"][condition] for condition in conditions},
        reaction_bounds=_bounds(context, conditions),
        joint_lambda=float(context["objective"]["joint_lambda"]),
        max_seconds=args.max_seconds,
        mip_gap=args.mip_gap,
        threads=args.threads,
        artifact_prefix=artifact_prefix,
        log_file=artifact_prefix.with_suffix(".gurobi.log"),
    )
    result_dict = result.to_dict()
    payload = {
        "status": result.status,
        "scientific_success": result.scientific_success,
        "schema_version": "metabolic_joint_checkpoint.v1",
        "context_sha256": context_sha,
        "study_accession": context["study_accession"],
        "sample_count": len(conditions),
        "solver": "gurobi",
        "joint_lambda": context["objective"]["joint_lambda"],
        "result": {
            "solver": "gurobi",
            "joint_lambda": context["objective"]["joint_lambda"],
            "active_tolerance": 1e-7,
            "conditions": conditions,
            "joint": result_dict["summaries"],
            "joint_active_union": result_dict["active_union"],
            "joint_active_union_size": result_dict["active_union_size"],
        },
        "instrumentation": result_dict,
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "claim_limit": (
            "Canonical result only when status=completed; partial incumbent is optimization "
            "diagnostic evidence and not a biological result."
        ),
    }
    if result.scientific_success:
        _atomic_write(args.output, payload)
        output = args.output
        exit_code = 0
    else:
        output = _attempt_path(args.attempt_dir, "joint")
        _atomic_write(output, payload)
        exit_code = 2
    print(json.dumps({"status": result.status, "output": str(output)}))
    return exit_code


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser()
    commands = result.add_subparsers(dest="command", required=True)
    for name, function in (("independent", independent), ("joint", joint)):
        command = commands.add_parser(name)
        command.add_argument("--context", type=Path, required=True)
        command.add_argument("--human-gem", type=Path, required=True)
        command.add_argument("--attempt-dir", type=Path, required=True)
        command.add_argument("--artifact-dir", type=Path, required=True)
        command.add_argument("--max-seconds", type=int, required=True)
        command.add_argument("--mip-gap", type=float, default=1e-4)
        command.add_argument("--threads", type=int, default=8)
        command.set_defaults(function=function)
    commands.choices["independent"].add_argument("--array-index", type=int, required=True)
    commands.choices["independent"].add_argument("--output-dir", type=Path, required=True)
    commands.choices["joint"].add_argument("--output", type=Path, required=True)
    return result


def main() -> int:
    args = parser().parse_args()
    if args.max_seconds <= 0 or args.threads <= 0 or not (0 <= args.mip_gap < 1):
        raise ValueError("invalid max-seconds, threads, or mip-gap")
    return int(args.function(args))


if __name__ == "__main__":
    raise SystemExit(main())

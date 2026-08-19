#!/usr/bin/env python3
"""Checkpointed recovery for cohort metabolic independent/joint CORNETO fits."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from dataclasses import asdict
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


def _existing_matches(path: Path, expected: dict[str, Any]) -> bool:
    if not path.is_file():
        return False
    value = _read_json(path)
    if value.get("status") != "completed" and value.get("status") != "prepared":
        raise ValueError(f"existing checkpoint is not complete: {path}")
    for key, wanted in expected.items():
        if value.get(key) != wanted:
            raise ValueError(
                f"existing checkpoint mismatch for {key}: {value.get(key)!r} != {wanted!r}"
            )
    return True


def _context(path: Path) -> tuple[dict[str, Any], str]:
    value = _read_json(path)
    if value.get("status") != "prepared":
        raise ValueError("scientific context is not prepared")
    return value, _sha256(path)


def _bounds(
    value: dict[str, Any], conditions: list[str]
) -> dict[str, dict[str, tuple[float, float]]]:
    result: dict[str, dict[str, tuple[float, float]]] = {}
    raw = value["reaction_bounds"]
    for condition in conditions:
        result[condition] = {
            reaction: (float(interval[0]), float(interval[1]))
            for reaction, interval in raw[condition].items()
        }
    return result


def prepare(args: argparse.Namespace) -> int:
    if _existing_matches(
        args.output,
        {"schema_version": "metabolic_checkpoint_context.v1", "study_accession": args.study},
    ):
        print(json.dumps({"status": "existing_valid", "output": str(args.output)}))
        return 0

    import cobra
    from run_corneto_14568_pilot import (
        _candidate_sets,
        _reaction_bounds,
        _read_expression,
        _read_manifest,
        _solver_choice,
    )

    rows = _read_manifest(args.manifest, study=args.study, primary_only=True)
    rows = rows[: args.expected_samples]
    if len(rows) != args.expected_samples:
        raise ValueError(f"expected {args.expected_samples} primary samples, found {len(rows)}")
    run_ids = [row["run_accession"] for row in rows]
    expression = _read_expression(args.expression, run_ids, "log1p_tpm")
    model = cobra.io.read_sbml_model(str(args.human_gem))
    solver, available_solvers, fallback_reason = _solver_choice("gurobi")
    if solver != "gurobi" or fallback_reason is not None:
        raise RuntimeError("checkpoint preparation requires explicit Gurobi without fallback")
    per_sample, audits, selected, selected_scores = _candidate_sets(
        model, expression, run_ids, max_candidates=25
    )
    reaction_bounds, growth_optima, growth_targets, bounds_skipped = _reaction_bounds(
        model, per_sample, run_ids, selected, 0.9
    )
    payload = {
        "status": "prepared",
        "schema_version": "metabolic_checkpoint_context.v1",
        "analysis": "cohort_metabolic_checkpoint_recovery",
        "study_accession": args.study,
        "primary_only": True,
        "sample_count": len(run_ids),
        "samples": rows,
        "conditions": run_ids,
        "expression": {
            "path": str(args.expression),
            "sha256": _sha256(args.expression),
            "transform": "log1p_tpm",
            "gene_count": int(expression.shape[0]),
        },
        "model": {
            "path": str(args.human_gem),
            "sha256": _sha256(args.human_gem),
            "reactions": len(model.reactions),
            "genes": len(model.genes),
            "biomass_id": "biomass_human",
        },
        "solver": {
            "requested": "gurobi",
            "used": solver,
            "available": available_solvers,
            "fallback_reason": None,
            "solve_fallback": None,
        },
        "candidate_selection": {
            "policy": "positive candidates ranked by median proposed upper bound; lowest retained",
            "max_candidates": 25,
            "selected_count": len(selected),
            "selected_reaction_ids": selected,
            "selected_median_proposed_upper": selected_scores,
            "audits_by_sample": audits,
            "bounds_clamped_to_human_gem": True,
            "bounds_skipped_as_empty_intersection": bounds_skipped,
        },
        "objective": {
            "biomass_coefficient": -1.0,
            "growth_fraction": 0.9,
            "growth_optima": growth_optima,
            "growth_targets": growth_targets,
            "independent_lambda": 0.1,
            "joint_lambda": 1.0,
        },
        "objectives": {run_id: {"biomass_human": -1.0} for run_id in run_ids},
        "reaction_bounds": reaction_bounds,
        "repo_commit": os.environ.get("HGSOC_CORNETO_REPO_COMMIT"),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "claim_limit": "Frozen response-blind execution context; no solved flux result.",
    }
    _atomic_write(args.output, payload)
    print(json.dumps({"status": "prepared", "samples": len(run_ids), "output": str(args.output)}))
    return 0


def solve_independent(args: argparse.Namespace) -> int:
    context, context_sha = _context(args.context)
    conditions = list(context["conditions"])
    if args.array_index < 0 or args.array_index >= len(conditions):
        raise IndexError(f"array index {args.array_index} outside 0..{len(conditions) - 1}")
    condition = conditions[args.array_index]
    output = args.output_dir / f"{args.array_index:03d}_{condition}.json"
    expected = {
        "schema_version": "metabolic_independent_checkpoint.v1",
        "condition": condition,
        "context_sha256": context_sha,
    }
    if _existing_matches(output, expected):
        print(json.dumps({"status": "existing_valid", "output": str(output)}))
        return 0

    import cobra

    from hgsoc_corneto.metabolic.joint_fba import solve_independent_sparse_fba

    model = cobra.io.read_sbml_model(str(args.human_gem))
    if _sha256(args.human_gem) != context["model"]["sha256"]:
        raise ValueError("Human-GEM SHA256 differs from the frozen context")
    summary = solve_independent_sparse_fba(
        model,
        condition=condition,
        objective=context["objectives"][condition],
        reaction_bounds=_bounds(context, [condition])[condition],
        independent_lambda=float(context["objective"]["independent_lambda"]),
        solver="gurobi",
    )
    payload = {
        "status": "completed",
        **expected,
        "study_accession": context["study_accession"],
        "array_index": args.array_index,
        "solver": "gurobi",
        "independent_lambda": context["objective"]["independent_lambda"],
        "solution": asdict(summary),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
    }
    _atomic_write(output, payload)
    print(json.dumps({"status": "completed", "condition": condition, "output": str(output)}))
    return 0


def solve_joint(args: argparse.Namespace) -> int:
    context, context_sha = _context(args.context)
    expected = {
        "schema_version": "metabolic_joint_checkpoint.v1",
        "context_sha256": context_sha,
    }
    if _existing_matches(args.output, expected):
        print(json.dumps({"status": "existing_valid", "output": str(args.output)}))
        return 0

    import cobra

    from hgsoc_corneto.metabolic.joint_fba import solve_joint_sparse_fba

    model = cobra.io.read_sbml_model(str(args.human_gem))
    if _sha256(args.human_gem) != context["model"]["sha256"]:
        raise ValueError("Human-GEM SHA256 differs from the frozen context")
    conditions = list(context["conditions"])
    result = solve_joint_sparse_fba(
        model,
        objectives={condition: context["objectives"][condition] for condition in conditions},
        reaction_bounds=_bounds(context, conditions),
        joint_lambda=float(context["objective"]["joint_lambda"]),
        solver="gurobi",
    )
    payload = {
        "status": "completed",
        **expected,
        "study_accession": context["study_accession"],
        "sample_count": len(conditions),
        "solver": "gurobi",
        "joint_lambda": context["objective"]["joint_lambda"],
        "result": result.to_dict(),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
    }
    _atomic_write(args.output, payload)
    print(
        json.dumps({"status": "completed", "samples": len(conditions), "output": str(args.output)})
    )
    return 0


def assemble(args: argparse.Namespace) -> int:
    context, context_sha = _context(args.context)
    conditions = list(context["conditions"])
    independent: list[dict[str, Any]] = []
    independent_jobs: list[str | None] = []
    for index, condition in enumerate(conditions):
        path = args.independent_dir / f"{index:03d}_{condition}.json"
        receipt = _read_json(path)
        expected = {
            "status": "completed",
            "schema_version": "metabolic_independent_checkpoint.v1",
            "condition": condition,
            "context_sha256": context_sha,
        }
        if any(receipt.get(key) != value for key, value in expected.items()):
            raise ValueError(f"invalid independent receipt: {path}")
        solution = receipt.get("solution")
        if not isinstance(solution, dict) or str(solution.get("status", "")).casefold() not in {
            "optimal",
            "optimal_inaccurate",
        }:
            raise ValueError(f"independent solution is not optimal: {path}")
        independent.append(solution)
        independent_jobs.append(receipt.get("slurm_job_id"))

    joint_receipt = _read_json(args.joint_receipt)
    if (
        joint_receipt.get("status") != "completed"
        or joint_receipt.get("schema_version") != "metabolic_joint_checkpoint.v1"
        or joint_receipt.get("context_sha256") != context_sha
    ):
        raise ValueError("joint checkpoint receipt is invalid")
    joint_result = joint_receipt.get("result")
    if not isinstance(joint_result, dict) or list(joint_result.get("conditions", [])) != conditions:
        raise ValueError("joint result conditions differ from frozen context")
    joint = joint_result.get("joint")
    if not isinstance(joint, list) or len(joint) != len(conditions):
        raise ValueError("joint result has the wrong number of condition summaries")
    if any(
        str(item.get("status", "")).casefold() not in {"optimal", "optimal_inaccurate"}
        for item in joint
    ):
        raise ValueError("at least one joint condition is not optimal")

    independent_union = sorted(
        {reaction for item in independent for reaction in item["active_by_flux"]}
    )
    joint_union = sorted({reaction for item in joint for reaction in item["active_by_flux"]})
    corneto = {
        "solver": "gurobi",
        "independent_lambda": context["objective"]["independent_lambda"],
        "joint_lambda": context["objective"]["joint_lambda"],
        "active_tolerance": joint_result["active_tolerance"],
        "conditions": conditions,
        "independent": independent,
        "joint": joint,
        "independent_active_union": independent_union,
        "joint_active_union": joint_union,
        "independent_active_union_size": len(independent_union),
        "joint_active_union_size": len(joint_union),
    }
    payload = {
        "status": "completed",
        "study_accession": context["study_accession"],
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "repo_commit": context.get("repo_commit"),
        "primary_only": True,
        "sample_count": context["sample_count"],
        "samples": context["samples"],
        "expression": {
            key: context["expression"][key] for key in ("path", "transform", "gene_count")
        },
        "model": {
            key: context["model"][key] for key in ("path", "reactions", "genes", "biomass_id")
        },
        "solver": context["solver"],
        "candidate_selection": context["candidate_selection"],
        "objective": context["objective"],
        "corneto": corneto,
        "checkpoint_provenance": {
            "schema_version": "metabolic_checkpoint_assembly.v1",
            "context_sha256": context_sha,
            "independent_job_ids": independent_jobs,
            "joint_job_id": joint_receipt.get("slurm_job_id"),
            "assembly_job_id": os.environ.get("SLURM_JOB_ID"),
        },
    }
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite canonical receipt: {args.output}")
    _atomic_write(args.output, payload)
    print(
        json.dumps({"status": "completed", "samples": len(conditions), "output": str(args.output)})
    )
    return 0


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser()
    commands = result.add_subparsers(dest="command", required=True)
    prep = commands.add_parser("prepare")
    prep.add_argument("--expression", type=Path, required=True)
    prep.add_argument("--manifest", type=Path, required=True)
    prep.add_argument("--human-gem", type=Path, required=True)
    prep.add_argument("--study", required=True)
    prep.add_argument("--expected-samples", type=int, required=True)
    prep.add_argument("--output", type=Path, required=True)
    prep.set_defaults(function=prepare)

    independent = commands.add_parser("independent")
    independent.add_argument("--context", type=Path, required=True)
    independent.add_argument("--human-gem", type=Path, required=True)
    independent.add_argument("--array-index", type=int, required=True)
    independent.add_argument("--output-dir", type=Path, required=True)
    independent.set_defaults(function=solve_independent)

    joint = commands.add_parser("joint")
    joint.add_argument("--context", type=Path, required=True)
    joint.add_argument("--human-gem", type=Path, required=True)
    joint.add_argument("--output", type=Path, required=True)
    joint.set_defaults(function=solve_joint)

    assembly = commands.add_parser("assemble")
    assembly.add_argument("--context", type=Path, required=True)
    assembly.add_argument("--independent-dir", type=Path, required=True)
    assembly.add_argument("--joint-receipt", type=Path, required=True)
    assembly.add_argument("--output", type=Path, required=True)
    assembly.set_defaults(function=assemble)
    return result


def main() -> int:
    args = parser().parse_args()
    return int(args.function(args))


if __name__ == "__main__":
    raise SystemExit(main())

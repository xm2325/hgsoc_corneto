"""Small metabolic fixtures that make order dependence falsifiable."""

from __future__ import annotations

from typing import Any

from hgsoc_corneto.metabolic.sequential import (
    CandidateConstraint,
    apply_sequential_constraints,
)


def parallel_pathway_model() -> Any:
    """Return a two-route model with a shared uptake capacity of ten."""

    try:
        from cobra import Metabolite, Model, Reaction
    except ImportError as error:  # pragma: no cover - exercised on Roihu
        raise RuntimeError("Install the 'metabolic' extra to build the toy model") from error

    model = Model("parallel_pathway_order_test")
    substrate = Metabolite("substrate_c", compartment="c")
    precursor = Metabolite("precursor_c", compartment="c")

    uptake = Reaction("UPTAKE")
    uptake.bounds = (0.0, 10.0)
    uptake.add_metabolites({substrate: 1.0})

    route_a = Reaction("ROUTE_A")
    route_a.bounds = (0.0, 1000.0)
    route_a.gene_reaction_rule = "gene_a"
    route_a.add_metabolites({substrate: -1.0, precursor: 1.0})

    route_b = Reaction("ROUTE_B")
    route_b.bounds = (0.0, 1000.0)
    route_b.gene_reaction_rule = "gene_b"
    route_b.add_metabolites({substrate: -1.0, precursor: 1.0})

    biomass = Reaction("BIOMASS")
    biomass.bounds = (0.0, 1000.0)
    biomass.add_metabolites({precursor: -1.0})

    model.add_reactions([uptake, route_a, route_b, biomass])
    model.objective = biomass
    return model


def _candidate(reaction_id: str) -> CandidateConstraint:
    gene = "gene_a" if reaction_id == "ROUTE_A" else "gene_b"
    return CandidateConstraint(
        reaction_id=reaction_id,
        category="single_gene_forward",
        genes=(gene,),
        expression_bound=2.0,
        proposed_lower=0.0,
        proposed_upper=2.0,
        reversible=False,
    )


def toy_order_benchmark(semantics: str = "bounds_safe") -> dict[str, Any]:
    """Run both candidate orders and summarize retained sets and FVA ranges."""

    try:
        from cobra.flux_analysis import flux_variability_analysis
    except ImportError as error:  # pragma: no cover - exercised on Roihu
        raise RuntimeError("Install the 'metabolic' extra to run the benchmark") from error

    runs: dict[str, Any] = {}
    for order in (("ROUTE_A", "ROUTE_B"), ("ROUTE_B", "ROUTE_A")):
        model = parallel_pathway_model()
        result = apply_sequential_constraints(
            model,
            [_candidate(reaction_id) for reaction_id in order],
            biomass_id="BIOMASS",
            growth_threshold=6.0,
            semantics=semantics,  # type: ignore[arg-type]
        )
        fva = flux_variability_analysis(
            model, reaction_list=["ROUTE_A", "ROUTE_B"], fraction_of_optimum=1.0
        )
        key = "_then_".join(order)
        runs[key] = {
            **result.to_dict(),
            "fva_at_optimum": {
                reaction_id: {
                    "minimum": float(fva.loc[reaction_id, "minimum"]),
                    "maximum": float(fva.loc[reaction_id, "maximum"]),
                }
                for reaction_id in ("ROUTE_A", "ROUTE_B")
            },
        }

    retained_sets = [set(run["retained_reactions"]) for run in runs.values()]
    intersection = retained_sets[0] & retained_sets[1]
    union = retained_sets[0] | retained_sets[1]
    return {
        "model": "two interchangeable routes with shared uptake capacity 10",
        "candidate_upper_bound": 2.0,
        "growth_threshold": 6.0,
        "semantics": semantics,
        "runs": runs,
        "retained_set_jaccard": len(intersection) / len(union) if union else 1.0,
        "order_dependence_observed": retained_sets[0] != retained_sets[1],
        "interpretation": (
            "Each constraint is feasible alone; the second makes growth fall to four and is "
            "reopened. Reversing order therefore reverses the retained reaction."
        ),
    }

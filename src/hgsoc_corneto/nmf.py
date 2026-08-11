"""Consensus diagnostics for label-blind, multi-start non-negative factorisation."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.cluster.hierarchy import cophenet, cut_tree, linkage
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import squareform
from sklearn.decomposition import NMF
from sklearn.metrics import silhouette_samples, silhouette_score


@dataclass(frozen=True)
class NmfRun:
    run_index: int
    seed: int
    reconstruction_error: float
    fit: float
    iterations: int
    converged: bool
    labels: np.ndarray


@dataclass(frozen=True)
class ConsensusNmfResult:
    rank: int
    consensus: np.ndarray
    labels: np.ndarray
    state_names: tuple[str, ...]
    assignment_stability: np.ndarray
    silhouette_by_sample: np.ndarray
    cophenetic_correlation: float
    dispersion: float
    consensus_sharpness: float
    average_silhouette: float
    runs: tuple[NmfRun, ...]
    best_w: np.ndarray
    best_h: np.ndarray


def select_top_mad(
    matrix: np.ndarray,
    gene_ids: tuple[str, ...],
    top_genes: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Select genes by descending median absolute deviation with lexical tie-breaking."""

    values = np.asarray(matrix, dtype=float)
    if values.ndim != 2 or values.shape[0] != len(gene_ids):
        raise ValueError("Gene matrix and identifiers have incompatible shapes")
    if top_genes < 1:
        raise ValueError("top_genes must be positive")
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("NMF input must be finite and non-negative")
    medians = np.median(values, axis=1)
    mad = np.median(np.abs(values - medians[:, None]), axis=1)
    eligible = np.flatnonzero(mad > 0)
    if len(eligible) < top_genes:
        raise ValueError(f"Only {len(eligible)} genes have positive MAD; requested {top_genes}")
    order = sorted(eligible, key=lambda index: (-mad[index], gene_ids[index]))
    selected = np.asarray(order[:top_genes], dtype=int)
    return selected, mad


def _deterministic_consensus_labels(
    consensus: np.ndarray,
    sample_ids: tuple[str, ...],
    rank: int,
) -> tuple[np.ndarray, np.ndarray, float]:
    distance = 1.0 - consensus
    np.fill_diagonal(distance, 0.0)
    condensed = squareform(distance, checks=True)
    tree = linkage(condensed, method="average")
    raw_labels = cut_tree(tree, n_clusters=rank).reshape(-1)
    unique = sorted(set(int(value) for value in raw_labels))
    if len(unique) != rank:
        raise RuntimeError(f"Consensus tree produced {len(unique)} clusters at rank {rank}")
    ordered = sorted(
        unique,
        key=lambda cluster: tuple(
            sorted(sample_ids[index] for index in np.flatnonzero(raw_labels == cluster))
        ),
    )
    label_map = {old: new for new, old in enumerate(ordered)}
    labels = np.asarray([label_map[int(value)] for value in raw_labels], dtype=int)
    coph = float(cophenet(tree, condensed)[0])
    return labels, distance, coph


def _align_labels(reference: np.ndarray, candidate: np.ndarray, rank: int) -> np.ndarray:
    overlap = np.zeros((rank, rank), dtype=int)
    for candidate_label in range(rank):
        for reference_label in range(rank):
            overlap[candidate_label, reference_label] = int(
                np.sum((candidate == candidate_label) & (reference == reference_label))
            )
    rows, columns = linear_sum_assignment(-overlap)
    mapping = {int(row): int(column) for row, column in zip(rows, columns, strict=True)}
    return np.asarray([mapping[int(value)] for value in candidate], dtype=int)


def run_consensus_nmf(
    matrix: np.ndarray,
    sample_ids: tuple[str, ...],
    *,
    rank: int,
    runs: int,
    seed_base: int,
    max_iter: int,
    tolerance: float,
) -> ConsensusNmfResult:
    """Fit random-start NMF models and summarize co-clustering stability."""

    values = np.asarray(matrix, dtype=float)
    if values.ndim != 2 or values.shape[1] != len(sample_ids):
        raise ValueError("NMF matrix must have genes as rows and samples as columns")
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("NMF input must be finite and non-negative")
    if rank < 2 or rank >= len(sample_ids):
        raise ValueError("rank must be at least 2 and smaller than the sample count")
    if runs < 2 or max_iter < 1 or tolerance <= 0:
        raise ValueError("Invalid NMF repetition or convergence parameters")

    denominator = float(np.square(values).sum())
    if denominator == 0:
        raise ValueError("NMF input is all zero")
    consensus_sum = np.zeros((len(sample_ids), len(sample_ids)), dtype=float)
    fitted_runs: list[NmfRun] = []
    best_error = float("inf")
    best_w: np.ndarray | None = None
    best_h: np.ndarray | None = None
    best_labels: np.ndarray | None = None
    for run_index in range(runs):
        seed = seed_base + rank * 100_000 + run_index
        model = NMF(
            n_components=rank,
            init="random",
            solver="cd",
            beta_loss="frobenius",
            random_state=seed,
            max_iter=max_iter,
            tol=tolerance,
            shuffle=False,
        )
        w = model.fit_transform(values)
        h = model.components_
        labels = np.argmax(h, axis=0).astype(int)
        error = float(model.reconstruction_err_)
        fit = 1.0 - error**2 / denominator
        fitted_runs.append(
            NmfRun(
                run_index=run_index,
                seed=seed,
                reconstruction_error=error,
                fit=fit,
                iterations=int(model.n_iter_),
                converged=bool(model.n_iter_ < max_iter),
                labels=labels,
            )
        )
        consensus_sum += labels[:, None] == labels[None, :]
        if error < best_error:
            best_error = error
            best_w = w.copy()
            best_h = h.copy()
            best_labels = labels.copy()

    consensus = consensus_sum / runs
    labels, distance, coph = _deterministic_consensus_labels(consensus, sample_ids, rank)
    off_diagonal = ~np.eye(len(sample_ids), dtype=bool)
    dispersion = float(4.0 * np.mean(np.square(consensus - 0.5)))
    sharpness = float(4.0 * np.mean(np.square(consensus[off_diagonal] - 0.5)))
    silhouette = silhouette_samples(distance, labels, metric="precomputed")
    average_silhouette = float(silhouette_score(distance, labels, metric="precomputed"))

    aligned_runs = tuple(
        NmfRun(
            run_index=run.run_index,
            seed=run.seed,
            reconstruction_error=run.reconstruction_error,
            fit=run.fit,
            iterations=run.iterations,
            converged=run.converged,
            labels=_align_labels(labels, run.labels, rank),
        )
        for run in fitted_runs
    )
    stability = np.mean(
        np.stack([run.labels == labels for run in aligned_runs], axis=0), axis=0
    )

    if best_w is None or best_h is None or best_labels is None:
        raise RuntimeError("No NMF solution was fitted")
    aligned_best = _align_labels(labels, best_labels, rank)
    component_to_state: dict[int, int] = {}
    for component in range(rank):
        matching_states = aligned_best[best_labels == component]
        if len(matching_states):
            component_to_state[component] = int(
                np.bincount(matching_states, minlength=rank).argmax()
            )
    if len(component_to_state) != rank or len(set(component_to_state.values())) != rank:
        overlap = np.zeros((rank, rank), dtype=int)
        for component in range(rank):
            for state in range(rank):
                overlap[component, state] = int(
                    np.sum((best_labels == component) & (labels == state))
                )
        rows, columns = linear_sum_assignment(-overlap)
        component_to_state = {
            int(component): int(state)
            for component, state in zip(rows, columns, strict=True)
        }
    state_to_component = {state: component for component, state in component_to_state.items()}
    component_order = [state_to_component[state] for state in range(rank)]
    state_names = tuple(f"tumour_state_{index + 1}" for index in range(rank))

    return ConsensusNmfResult(
        rank=rank,
        consensus=consensus,
        labels=labels,
        state_names=state_names,
        assignment_stability=stability,
        silhouette_by_sample=np.asarray(silhouette, dtype=float),
        cophenetic_correlation=coph,
        dispersion=dispersion,
        consensus_sharpness=sharpness,
        average_silhouette=average_silhouette,
        runs=aligned_runs,
        best_w=best_w[:, component_order],
        best_h=best_h[component_order, :],
    )

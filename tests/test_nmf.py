import numpy as np

from hgsoc_corneto.nmf import run_consensus_nmf, select_top_mad


def test_top_mad_selection_has_deterministic_ties() -> None:
    matrix = np.asarray(
        [
            [0, 2, 0, 2],
            [0, 1, 0, 1],
            [0, 2, 0, 2],
            [1, 1, 1, 1],
        ],
        dtype=float,
    )
    selected, mad = select_top_mad(matrix, ("gene_b", "gene_c", "gene_a", "gene_d"), 2)
    assert selected.tolist() == [2, 0]
    assert mad[selected].tolist() == [1.0, 1.0]


def test_consensus_nmf_recovers_two_separated_blocks() -> None:
    matrix = np.asarray(
        [
            [12, 11, 10, 0, 0, 0],
            [9, 10, 11, 0, 0, 0],
            [7, 8, 9, 0, 0, 0],
            [0, 0, 0, 10, 11, 12],
            [0, 0, 0, 11, 10, 9],
            [0, 0, 0, 8, 9, 7],
        ],
        dtype=float,
    )
    result = run_consensus_nmf(
        matrix,
        ("S1", "S2", "S3", "S4", "S5", "S6"),
        rank=2,
        runs=12,
        seed_base=100,
        max_iter=1000,
        tolerance=1e-5,
    )
    assert result.state_names == ("tumour_state_1", "tumour_state_2")
    assert result.labels.tolist() == [0, 0, 0, 1, 1, 1]
    assert result.consensus_sharpness > 0.99
    assert result.average_silhouette > 0.99
    assert result.best_w.shape == (6, 2)
    assert result.best_h.shape == (2, 6)

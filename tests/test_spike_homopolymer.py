"""Unit tests for the produceBCmat3 ports: CC-level spike detection +
homopolymer filtering.

Verifies:
  - identify_spike_ccs flags whole CC as spike when any member ends with tag
  - has_homopolymer_run / filter_homopolymers behave like findhomopolymers.m
  - split_spike_ins_by_cc uses the CC flag (not representative check)
"""

import sys
from pathlib import Path

PROJ_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJ_ROOT))

import numpy as np

from mapseq2.errorcorrect import identify_spike_ccs
from mapseq2.matrix import (
    filter_homopolymers,
    has_homopolymer_run,
    split_spike_ins,
    split_spike_ins_by_cc,
)


def test_homopolymer_run():
    assert not has_homopolymer_run("ACGTACGT", 7)
    assert not has_homopolymer_run("AAAAAA" + "CGT", 7)
    assert has_homopolymer_run("AAAAAAA" + "CGT", 7)
    assert has_homopolymer_run("CGT" + "TTTTTTTT", 7)
    assert has_homopolymer_run("CG" + "CCCCCCCC" + "GT", 7)
    assert not has_homopolymer_run("CGT", 7)
    assert not has_homopolymer_run("", 7)
    assert not has_homopolymer_run("AAAAAAAA", 0)
    assert not has_homopolymer_run("AAAAAAAA", 1)
    print("homopolymer_run tests passed")


def test_identify_spike_ccs():
    pooled = [
        "AAAAAAAAAAAAAAAAAAAAAAAAACGTCAGTC",  # spike (last 8 = CGTCAGTC)
        "AAAAAAAAAAAAAAAAAAAAAAAAACGTCAGTA",  # 1bp variant of spike
        "AAAAAAAAAAAAAAAAAAAAAAAAACGTCAATC",  # another spike variant
        "GGGGTTTTGGGGTTTTGGGGTTTTGGGGTT",     # unrelated barcode
    ]
    spike_tag = "CGTCAGTC"

    # CC 0 contains pristine spike + 2 variants (real-world clustering scenario)
    # CC 1 contains the unrelated barcode alone
    global_labels = np.array([0, 0, 0, 1], dtype=np.int64)
    is_spike = identify_spike_ccs(pooled, global_labels, spike_tag)
    assert is_spike.tolist() == [True, False]
    print("CC 0 = spike (one member matches tag), CC 1 = real: ✓")

    # If pristine spike is missing (only variants in CC), it would NOT be flagged
    pooled_no_pristine = pooled[1:]
    labels_no_pristine = np.array([0, 0, 1], dtype=np.int64)
    is_spike_v = identify_spike_ccs(pooled_no_pristine, labels_no_pristine, spike_tag)
    assert is_spike_v.tolist() == [False, False]
    print("when CC has only variants (no pristine tag), not flagged (correct)")


def test_split_by_cc_vs_representative():
    """The whole reason CC-based detection exists: when CC representative is a
    variant (not the pristine tag), rep-based check misses it but CC-check
    catches it."""
    matrix = np.array([
        [100, 50, 80],
        [10, 20, 30],
    ], dtype=np.int64)
    rep_pristine = ["AAAAAAAAAAAAAAAAAAAAAAAAACGTCAGTC", "GGGGTTTTGGGGTTTTGGGGTTTTGGGGTT"]
    rep_variant = ["AAAAAAAAAAAAAAAAAAAAAAAAACGTCAGTA", "GGGGTTTTGGGGTTTTGGGGTTTTGGGGTT"]

    r1, _, s1, _ = split_spike_ins(matrix, rep_variant, "CGTCAGTC")
    assert r1.shape[0] == 2 and s1.shape[0] == 0
    print("rep-based check: variant rep → NOT classified as spike (would lose accuracy)")

    is_spike_cc = np.array([True, False])
    r2, _, s2, _ = split_spike_ins_by_cc(matrix, rep_variant, is_spike_cc)
    assert r2.shape[0] == 1 and s2.shape[0] == 1
    assert np.array_equal(s2[0], matrix[0])
    print("CC-based check: variant rep flagged via CC → correctly classified")


def test_filter_homopolymers():
    matrix = np.array([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
    ], dtype=np.int64)
    bcs = [
        "ACGTACGTACGTACGTACGTACGTACGTACGT",
        "AAAAAAACGTACGTACGTACGTACGTACGTAC",  # 7 A's in a row
        "CGTACGTACGTACGTACGTACGTACGTACGTA",
    ]
    is_spike = np.array([False, True, False])

    new_mat, new_bcs, keep, (new_spike,) = filter_homopolymers(matrix, bcs, 7, extra=(is_spike,))
    assert new_mat.shape == (2, 3)
    assert len(new_bcs) == 2
    assert keep.tolist() == [True, False, True]
    assert new_spike.tolist() == [False, False]
    print("homopolymer filter drops row with 7-A run; extra arrays stay aligned")

    new_mat0, _, keep0, _ = filter_homopolymers(matrix, bcs, 0)
    assert new_mat0.shape == matrix.shape
    print("min_run=0 disables filter")


def main():
    test_homopolymer_run()
    test_identify_spike_ccs()
    test_split_by_cc_vs_representative()
    test_filter_homopolymers()
    print("\nALL produceBCmat3 PORT TESTS PASSED")


if __name__ == "__main__":
    main()

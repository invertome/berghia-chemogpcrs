"""Tests for the SLURM array-range planner (bead fxx).

Stages 04/05 hardcode `#SBATCH --array=0-999%50`, which cannot cover a
full-557 manifest. The submit wrapper must size the array from OG_COUNT
(and detect the MaxArraySize ceiling) instead.
"""
import pytest

from array_plan import array_spec, chunk_ranges, fits_single_array

MAX_ARRAY_SIZE = 10001   # Unity: max array size 10001 => max index 10000


def test_spec_covers_every_task_index():
    # n tasks -> indices 0..n-1
    assert array_spec(1, throttle=50) == "0-0%50"
    assert array_spec(1000, throttle=50) == "0-999%50"
    assert array_spec(1743, throttle=20) == "0-1742%20"


def test_spec_without_throttle_has_no_suffix():
    assert array_spec(500) == "0-499"


def test_spec_rejects_non_positive_task_count():
    # A zero-task manifest must fail loudly, never emit a bogus range.
    with pytest.raises(ValueError):
        array_spec(0, throttle=50)


def test_exactly_max_array_size_still_fits_one_array():
    # size 10001 => indices 0..10000, the highest legal index
    assert fits_single_array(MAX_ARRAY_SIZE, MAX_ARRAY_SIZE) is True
    assert array_spec(MAX_ARRAY_SIZE, throttle=50) == "0-10000%50"


def test_one_over_max_array_size_needs_chunking():
    assert fits_single_array(MAX_ARRAY_SIZE + 1, MAX_ARRAY_SIZE) is False


def test_chunk_ranges_are_contiguous_and_cover_all_tasks():
    ranges = chunk_ranges(MAX_ARRAY_SIZE + 1, MAX_ARRAY_SIZE)
    assert ranges == [(0, 10000), (10001, 10001)]

    ranges = chunk_ranges(25000, MAX_ARRAY_SIZE)
    assert ranges == [(0, 10000), (10001, 20001), (20002, 24999)]
    # contiguous, no gaps/overlaps, covers 0..n-1
    assert ranges[0][0] == 0
    assert ranges[-1][1] == 24999
    for (_, prev_end), (next_start, _) in zip(ranges, ranges[1:]):
        assert next_start == prev_end + 1
    # every chunk is within the ceiling
    assert all((end - start + 1) <= MAX_ARRAY_SIZE for start, end in ranges)


def test_chunk_ranges_single_chunk_when_it_fits():
    assert chunk_ranges(1000, MAX_ARRAY_SIZE) == [(0, 999)]

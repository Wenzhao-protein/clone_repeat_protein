from __future__ import annotations

import importlib.util
import random
from pathlib import Path

import numpy as np
import pandas as pd

from hurdler.constants import AMINO_ACIDS, PLASMIDS


def _load_script():
    path = Path(__file__).parents[1] / "scripts" / "run_success_landscape_single_files.py"
    specification = importlib.util.spec_from_file_location("single_file_success_landscape", path)
    assert specification is not None and specification.loader is not None
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def test_short_task_partition_is_complete_and_nonoverlapping():
    module = _load_script()
    for maximum in range(1, 6):
        tasks = [task for task in module._short_tasks(maximum) if task[0] == maximum]
        assert sum(20 ** (maximum - len(prefix)) for _, prefix in tasks) == 20**maximum
        assert len({prefix for _, prefix in tasks}) == len(tasks)
    motifs = [
        motif
        for length, prefix in module._short_tasks(3)
        if length == 3
        for motif in module._motifs(length, prefix)
    ]
    assert len(motifs) == 8000
    assert len(set(motifs)) == 8000


def test_random_task_generation_preserves_frozen_seed_blocks():
    module = _load_script()
    tasks = list(module._random_tasks(6, 7, 3))
    assert len(tasks) == 2 * len(PLASMIDS)

    six_generator = random.Random(module.SIX_AA_SEED)
    expected_six = [
        "".join(six_generator.choice(AMINO_ACIDS) for _ in range(6))
        for _ in range(3)
    ]
    assert tasks[0] == (6, PLASMIDS[0], module.SIX_AA_SEED, expected_six)

    seven_generator = random.Random(module.LEGACY_SEED)
    expected_seven = [
        "".join(seven_generator.choice(AMINO_ACIDS) for _ in range(7))
        for _ in range(3)
    ]
    assert tasks[len(PLASMIDS)] == (7, PLASMIDS[0], module.LEGACY_SEED, expected_seven)


def test_output_schemas_have_one_column_per_maintained_plasmid():
    module = _load_script()
    schema_names = set(module._short_schema().names)
    for plasmid in PLASMIDS:
        assert f"{plasmid}_success" in schema_names
    assert {"module_length", "plasmid", "module", "success"}.issubset(
        module.RANDOM_SCHEMA.names
    )


def test_three_copy_scan_is_the_active_cli_default(monkeypatch):
    module = _load_script()
    monkeypatch.setattr(module.sys, "argv", ["run_success_landscape_single_files.py"])
    arguments = module.parse_args()
    assert arguments.scan_copies == 3
    assert "run06_three_copy_16core" in str(arguments.short_output)
    assert arguments.figure_dir.name == "scan_3x"


def test_success_plot_uses_requested_near_square_style(tmp_path, monkeypatch):
    module = _load_script()
    frame = pd.DataFrame(
        [
            {"module_length": length, "plasmid": plasmid, "success_rate": 0.5}
            for plasmid in PLASMIDS
            for length in (1, 50, 60)
        ]
    )
    monkeypatch.setattr(module.matplotlib.figure.Figure, "savefig", lambda *a, **k: None)
    monkeypatch.setattr(module.plt, "close", lambda *a, **k: None)
    module.plot_success_curve(frame, tmp_path, file_stem="style_test")
    figure = module.plt.gcf()
    axis = figure.axes[0]
    assert tuple(figure.get_size_inches()) == (6.0, 5.0)
    assert axis.get_title() == "3-mer Probability vs Sequence Length"
    assert axis.get_xlabel() == "Sequence Length"
    assert axis.get_ylabel() == "Probability (%)"
    assert tuple(axis.get_xlim()) == (1.0, 50.0)
    assert all(max(line.get_xdata()) == 50 for line in axis.lines)
    assert [line.get_label() for line in axis.lines] == [
        "pET-28a(+)",
        "pET-28a(+)_start_codon",
        "pGEX-4T-1",
        "pMAL-c5X",
        "pUC18",
        "pQE-3",
        "pCold_I",
        "pET-21a(+)",
    ]


def test_historical_matcher_uses_doubled_module_for_cross_boundary_hit():
    module = _load_script()
    masks = np.zeros(8000 * 8000, dtype=np.uint8)
    left = module._encode_three_mer("FGA")
    right = module._encode_three_mer("EFG")
    masks[left * 8000 + right] = 1
    present = np.zeros(8000, dtype=np.bool_)
    present[[left, right]] = True
    index = module.HistoricalPatternIndex(masks, present, {})

    # FGA starts at position 4 and EFG at position 9 only after ACDEFG is
    # doubled. Their start distance is 5, which is valid for a 6AA module.
    assert module.historical_success_mask("ACDEFG", index) == 1


def test_historical_pattern_population_matches_committed_notebook_inputs():
    module = _load_script()
    index = module.build_historical_pattern_index(Path(__file__).parents[1])
    expected = {
        "pGEX-4T-1": 895_481,
        "pMAL-c5X": 736_353,
        "pET-21a(+)": 718_033,
        "pET-28a(+)": 545_149,
        "pET-28a(+)_start_codon": 723_332,
        "pCold_I": 1_021_688,
        "pUC18": 1_292_709,
        "pQE-3": 762_145,
    }
    for bit, plasmid in enumerate(PLASMIDS):
        assert int(np.count_nonzero(index.masks & (1 << bit))) == expected[plasmid]


def test_three_copy_scan_never_loses_two_copy_matches_and_can_add_boundary_hit():
    module = _load_script()
    index = module.build_historical_pattern_index(Path(__file__).parents[1])
    generator = random.Random(20260810)
    for length in (1, 2, 3, 5, 6, 7, 15, 31, 60):
        for _ in range(20):
            sequence = "".join(generator.choice(AMINO_ACIDS) for _ in range(length))
            two_copy = module.historical_success_mask(sequence, index, scan_copies=2)
            three_copy = module.historical_success_mask(sequence, index, scan_copies=3)
            assert two_copy & ~three_copy == 0

    # This 6AA module only becomes compatible when the third copy completes
    # the cyclic window that starts at the final residue of the first copy.
    assert module.historical_success_mask("VRFGEP", index, scan_copies=2) == 0
    assert module.historical_success_mask("VRFGEP", index, scan_copies=3) != 0

from hap_counter.plotting import plot_discrepancy_histogram


def test_plot_discrepancy_histogram_writes_nonempty_file(tmp_path):
    output_path = tmp_path / "discrepancy.png"
    fractions_by_group = {
        ("chr1", 1): [0.0, 0.1, 0.2, 0.05],
        ("chr1", 2): [0.3, 0.4],
        ("chr2", 1): [0.0, 0.0, 0.5],
    }

    plot_discrepancy_histogram(fractions_by_group, "some summary text", str(output_path))

    assert output_path.exists()
    assert output_path.stat().st_size > 0


def test_plot_discrepancy_histogram_handles_empty_input(tmp_path):
    output_path = tmp_path / "empty.png"

    plot_discrepancy_histogram({}, "no data", str(output_path))

    assert output_path.exists()
    assert output_path.stat().st_size > 0

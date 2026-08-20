from hap_counter.cli import _sample_mismatch_warning


def test_sample_mismatch_warning_none_when_bam_has_no_sample_names():
    message = _sample_mismatch_warning([], "hg002")
    assert message is not None
    assert "no sample name" in message.lower()
    assert "hg002" in message


def test_sample_mismatch_warning_none_when_names_match():
    assert _sample_mismatch_warning(["hg002"], "hg002") is None


def test_sample_mismatch_warning_none_when_vcf_sample_among_multiple_bam_samples():
    assert _sample_mismatch_warning(["hg002", "other"], "hg002") is None


def test_sample_mismatch_warning_when_names_differ():
    message = _sample_mismatch_warning(["NA12878"], "hg002")
    assert message is not None
    assert "NA12878" in message
    assert "hg002" in message

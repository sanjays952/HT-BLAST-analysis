"""Tests for the hit_filter module."""

from src.blast_project.hit_filter import filter_hits


class MockHsp:
    """Mock HSP object for testing."""

    def __init__(self, expect, query_start, query_end, identities, align_length, score):
        self.expect = expect
        self.query_start = query_start
        self.query_end = query_end
        self.identities = identities
        self.align_length = align_length
        self.score = score


class MockAlignment:
    """Mock Alignment object for testing."""

    def __init__(self, title, accession, hsps):
        self.title = title
        self.accession = accession
        self.hsps = hsps


class MockBlastRecord:
    """Mock BLAST record for testing."""

    def __init__(self, alignments):
        self.alignments = alignments


def test_filter_hits_returns_top_three():
    """Should return at most 3 hits."""
    hsps = [
        MockHsp(expect=1e-10, query_start=1, query_end=80, identities=72, align_length=80, score=300)
    ]
    alignments = [
        MockAlignment(title="gi|123|ref|XP_001| Species A", accession="XP_001", hsps=hsps),
        MockAlignment(title="gi|124|ref|XP_002| Species B", accession="XP_002", hsps=hsps),
        MockAlignment(title="gi|125|ref|XP_003| Species C", accession="XP_003", hsps=hsps),
        MockAlignment(title="gi|126|ref|XP_004| Species D", accession="XP_004", hsps=hsps),
    ]
    record = MockBlastRecord(alignments)

    result = filter_hits(record, "A" * 100)
    assert len(result) <= 3


def test_filter_hits_no_matches():
    """Should return empty list when nothing passes filters."""
    hsps = [
        MockHsp(expect=10, query_start=1, query_end=5, identities=1, align_length=100, score=10)
    ]
    alignments = [
        MockAlignment(title="gi|123|ref|XP_001| Bad Hit", accession="XP_001", hsps=hsps),
    ]
    record = MockBlastRecord(alignments)

    result = filter_hits(record, "A" * 100)
    assert result == []

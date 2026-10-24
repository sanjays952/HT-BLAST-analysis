"""Filter and rank BLAST hits by E-value, coverage, and identity."""


def filter_hits(
    blast_record,
    query_sequence,
    evalue_threshold=1e-5,
    query_cover_range=(40, 100),
    identity_range=(70, 100),
    score_threshold=200,
):
    """Filter BLAST hits and return the top three by E-value."""
    filtered_hits = []

    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if (
                hsp.expect <= evalue_threshold
                and query_cover_range[0]
                <= (hsp.query_end - hsp.query_start + 1) / len(query_sequence) * 100
                <= query_cover_range[1]
                and identity_range[0]
                <= (hsp.identities / hsp.align_length) * 100
                <= identity_range[1]
                and hsp.score > score_threshold
            ):
                species_name = alignment.title.split("|")[-1].strip()
                percentage_identity = (hsp.identities / hsp.align_length) * 100
                query_cover = (
                    (hsp.query_end - hsp.query_start + 1) / len(query_sequence) * 100
                )

                filtered_hits.append(
                    {
                        "Species Name": species_name,
                        "Accession": alignment.accession,
                        "Score": hsp.score,
                        "E-Value": hsp.expect,
                        "Query Cover": query_cover,
                        "Percentage Identity": percentage_identity,
                    }
                )

    # Sort hits by E-value (ascending), query cover (descending),
    # and percentage identity (descending)
    sorted_hits = sorted(
        filtered_hits,
        key=lambda x: (x["E-Value"], -x["Query Cover"], -x["Percentage Identity"]),
    )

    return sorted_hits[:3]

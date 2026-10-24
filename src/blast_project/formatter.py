"""Format and write BLAST hit results."""

from prettytable import PrettyTable


def format_hits_table(hits):
    """Create a PrettyTable from a list of hit dicts."""
    table = PrettyTable()
    table.field_names = [
        "Species Name",
        "Accession",
        "Score",
        "E-Value",
        "Query Cover",
        "Percentage Identity",
    ]
    for hit in hits:
        table.add_row(
            [
                hit["Species Name"],
                hit["Accession"],
                hit["Score"],
                hit["E-Value"],
                f"{hit['Query Cover']:.2f}%",
                f"{hit['Percentage Identity']:.2f}%",
            ]
        )
    return table


def write_results(output_path, header, table=None):
    """Append results for a single query to the output file."""
    with open(output_path, "a") as f:
        if table is None:
            f.write(header)
            f.write("Has no result\n\n")
        else:
            f.write(str(table))
            f.write("\n")
            f.write(header)
            f.write("\n")

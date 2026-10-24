"""Entry point for BLAST sequence analysis."""

from src.blast_project.blast_runner import parse_fasta, run_blastx
from src.blast_project.formatter import format_hits_table, write_results
from src.blast_project.hit_filter import filter_hits

OUTPUT_FILE = "output/top_three_blast_hits_table.txt"


def main():
    queries = parse_fasta("data/input1.fasta")

    for header, sequence in queries:
        blast_record = run_blastx(sequence)
        top_hits = filter_hits(blast_record, sequence)

        if not top_hits:
            print(f"No matching result for: {header}")
            write_results(OUTPUT_FILE, header)
        else:
            table = format_hits_table(top_hits)
            print(table)
            write_results(OUTPUT_FILE, header, table)


if __name__ == "__main__":
    main()

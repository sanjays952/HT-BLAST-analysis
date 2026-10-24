"""Parse FASTA files and run NCBI BLASTx queries."""

from Bio.Blast import NCBIWWW, NCBIXML


def parse_fasta(filepath):
    """Parse a FASTA file into a list of (header, sequence) tuples."""
    queries = []
    current_header = None
    current_sequence = []

    with open(filepath) as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    queries.append((current_header, "".join(current_sequence)))
                current_header = line
                current_sequence = []
            else:
                current_sequence.append(line)

    if current_header is not None:
        queries.append((current_header, "".join(current_sequence)))

    return queries


def run_blastx(sequence, database="refseq_protein"):
    """Run a remote BLASTx query and return the parsed record."""
    result_handle = NCBIWWW.qblast("blastx", database, sequence)
    blast_record = NCBIXML.read(result_handle)
    result_handle.close()
    return blast_record

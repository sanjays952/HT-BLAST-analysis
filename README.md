# Blast Project

NCBI BLASTx sequence analysis tool that:
- Reads FASTA sequences from input files
- Runs remote BLASTx against the RefSeq Protein database
- Filters hits by E-value, query coverage, identity, and score
- Outputs top 3 hits per query as a formatted table

## Setup

```bash
python -m venv .venv
.venv\Scripts\activate
pip install -e .

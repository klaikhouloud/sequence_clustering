
A Python script for clustering and filtering DNA sequences using VSEARCH and biological criteria like sequence complexity and GC content.

this script can :

1. Clusters sequences from a FASTA file using [VSEARCH](https://github.com/torognes/vsearch)
2. Filters clusters based on:
   - Minimum copy number
   - Overlapping genomic positions
   - Sequence complexity (LCC)
   - GC content
   - Diversity of flanking regions
3. Outputs filtered clusters as a new FASTA file


# Requirements

- Python ≥ 3.6
- [Biopython](https://biopython.org/) (`pip install biopython`)
- [VSEARCH](https://github.com/torognes/vsearch) (binary must be available locally)

# File structure

You need to prepare the following files:

- `file_candidates_fasta`: FASTA file containing sequences to cluster
- `file_temp_cluster_dir`: Directory where VSEARCH cluster files will be stored
- `file_family_fasta`: Output file for filtered clusters
- `file_temp_cluster`: Prefix for VSEARCH cluster files

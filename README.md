# fasta2nex

fasta2nex is a small python script that converts nucleotide alignments from FASTA
format to codon-partitioned Nexus format.

## Requirements

fasta2nex must be run on Linux system. It requires the following programs:

* ncbi-blast+ (can be installed with `apt install ncbi-blast+`)
* EMBOSS (can be installed [here](http://emboss.open-bio.org/html/use/ch02s05.html))

fasta2nex also requires the following Python packages:

* biopython 
* subprocess
* re 
* argparse
* datetime

### Usage

fasta2nex takes the following as input:

1. A multi-FASTA file of a gene
2. An aligned FASTA file of a gene 
3. A blast database of a gene

The gene name must be the third item in the FASTA header, separated by spaces.
For example:

```fasta
>Eurycea_waterlooensis TRINITY_283746 shh
```

You can run this script using Python or make it executable like this:

```bash
chmod +x fasta2nex.py
```

The command follows this general structure. The path to the multi-fasta file goes
after `-m`; the path to the alignment goes after `-a`, anf the path to the blast
database goes after `-d`.

If not executable:

```bash
python3 fasta2nex -m path/to/multi.fasta -a path/to/aln.fasta -d path/to/blastdb
```

If executable:

```bash
./fasta2nex -m path/to/multi.fasta -a path/to/aln.fasta -d path/to/blastdb
```

Output will be written to a file called `gene_name.nex`.
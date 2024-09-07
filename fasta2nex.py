#!/bin/env python

from Bio import Seq, SeqIO, SeqRecord
import subprocess as subp
import re
import argparse
from datetime import datetime

print("\n\
███████╗░█████╗░░██████╗████████╗░█████╗░██████╗░███╗░░██╗███████╗██╗░░██╗\n\
██╔════╝██╔══██╗██╔════╝╚══██╔══╝██╔══██╗╚════██╗████╗░██║██╔════╝╚██╗██╔╝\n\
█████╗░░███████║╚█████╗░░░░██║░░░███████║░░███╔═╝██╔██╗██║█████╗░░░╚███╔╝░\n\
██╔══╝░░██╔══██║░╚═══██╗░░░██║░░░██╔══██║██╔══╝░░██║╚████║██╔══╝░░░██╔██╗░\n\
██║░░░░░██║░░██║██████╔╝░░░██║░░░██║░░██║███████╗██║░╚███║███████╗██╔╝╚██╗\n\
╚═╝░░░░░╚═╝░░╚═╝╚═════╝░░░░╚═╝░░░╚═╝░░╚═╝╚══════╝╚═╝░░╚══╝╚══════╝╚═╝░░╚═╝\n")

print("FASTA2NEX: A small python utility to convert FASTA files to codon partition files.\n")
print("Requirements: Python 3, local blastp installation, multi fasta file and alignment, blast database.\n")
print("Example: ./fasta2nex.py -m multifasta -a alignment -d blast_database\n")

parser = argparse.ArgumentParser()

parser.add_argument("-m", "--multi", help = "Path to multifasta file.")
parser.add_argument("-a", "--aln", help = "Path to alignment file.")
parser.add_argument("-d", "--database", help = "Path to BLAST database.")

args = parser.parse_args()

if args.multi:
	print(f"Path to multi FASTA file: {args.multi}")
else:
	print("ERROR: paths to files/databases nonexistent. Exiting.\n")
	exit()
if args.aln:
	print(f"Path to alignment file: {args.aln}")
else:
        print("ERROR: paths to files/databases nonexistent. Exiting.\n")
        exit()
if args.database:
	print(f"Path to BLAST database: {args.database}\n")
else:
        print("ERROR: paths to files/databases nonexistent. Exiting.\n")
        exit()

current_time = datetime.now()
print("Begun at", current_time)

multifasta = SeqIO.parse(args.multi, "fasta")
longest_record = max(multifasta, key=len)
taxon_name = longest_record.id
gene_name = longest_record.description.split(" ")[2]

candidate_orfs = [longest_record[i:].translate(to_stop=True) for i in range(3, len(longest_record))]
top_orfs = round(len(candidate_orfs) * 0.1)

top_candidate_orfs = sorted(candidate_orfs, key=len)[-top_orfs:]

print("Conducting BLAST search. Please hang tight.\n")

iter = 0
evalue_winner = 10000000
len_winner = 0
for record in top_candidate_orfs:
	record.id = "Sequence" + str(iter)

	record_out = f".{gene_name}.{record.id}.fasta"
	SeqIO.write(record, record_out, "fasta")

	db = f"{args.database}"
	out = f".{gene_name}.{record.id}.outfmt6"
	subp.run(f"blastp -query {record_out} -db {db} -out {out} -outfmt 6 -evalue 1e-20",
		shell=True)

	with open(out, "r") as file:
		blast_record = file.read()

	if blast_record == "":
		iter += 1
		continue

	percent_id = blast_record.split("\t")[2]
	evalue = blast_record.split("\t")[10]

	if len(record.seq) >= len_winner and float(evalue) < evalue_winner:
		evalue_winner = float(evalue)
		len_winner = len(record.seq)
		overall_winner = record

	iter += 1

print("Best AA sequence is:", overall_winner.seq)
print("Length:", len_winner, "\n")

longest_record_pep1 = longest_record.seq.translate(to_stop=True)
longest_record_pep1.id = "frame1"

longest_record_pep2 = longest_record.seq[1:].translate(to_stop=True)
longest_record_pep2.id = "frame2"

longest_record_pep3 = longest_record.seq[2:].translate(to_stop=True)
longest_record_pep3.id = "frame3"

if overall_winner.seq in longest_record_pep1:
	frame_correction = 0
	print("AA sequence matches frame 1.")
elif overall_winner.seq in longest_record_pep2:
	frame_correction = 1
	print("AA sequence matches frame 2.")
elif overall_winner.seq in longest_record_pep3:
	frame_correction = 2
	print("AA sequence matches frame 3.")
else:
	print("ERROR: AA sequence doesn't match any frame. Exiting.\n")
	exit

aln = SeqIO.parse(args.aln, "fasta")
for record in aln:
	if taxon_name in record.id:
		aln_rec_of_int = record

gap = 0
for pos in aln_rec_of_int.seq:
	if pos == "-":
		gap += 1
	else:
		break

gap_correction = gap % 3

overall_correction = (frame_correction + gap_correction) % 3

begin_pos = 1 + overall_correction

print(gap, "leading gaps detected.")
print("Overall Correction ==", overall_correction, "\n")

end_pos1 = (len(aln_rec_of_int.seq) - (len(aln_rec_of_int.seq) % 3)) - 2
end_pos2 = (len(aln_rec_of_int.seq) - (len(aln_rec_of_int.seq) % 3)) - 1
end_pos3 = (len(aln_rec_of_int.seq) - (len(aln_rec_of_int.seq) % 3))

nexus_pos1 = "charset " + gene_name + "pos1" + " = " + str(begin_pos) + "-" + str(end_pos1) + "\\3;\n"
nexus_pos2 = "charset " + gene_name + "pos2" + " = " + str(begin_pos + 1) + "-" + str(end_pos2) + "\\3;\n"
nexus_pos3 = "charset " + gene_name + "pos3" + " = " + str(begin_pos + 2) + "-" + str(end_pos3) + "\\3;\n"

nexus_header = "\n\n#NEXUS\n\n\nbegin sets;\n"
nexus_footer = "end;"

bottom_nexus_file = nexus_header + nexus_pos1 + nexus_pos2 + nexus_pos3 + nexus_footer

nexus_bottom_name = f"{gene_name}_bottom.nex"
with open(nexus_bottom_name, "w") as file:
	file.write(bottom_nexus_file)

subp.run(f"seqret -sequence {args.aln} -outseq {args.aln}.nex -osformat nexus",
    shell = True)

subp.run(f"cat {args.aln}.nex {nexus_bottom_name} > {gene_name}.nex",
    shell = True)

print(f"Output written to {gene_name}.nex")

current_time = datetime.now()
print("Ended at", current_time)

subp.run(f"rm .*fasta .*outfmt6 {args.aln}.nex {nexus_bottom_name}",
    shell=True)

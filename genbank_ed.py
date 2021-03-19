#!/usr/bin/env python3
# coding: utf-8

from Bio import Entrez, SeqIO
from Bio.SeqUtils.CheckSum import seguid
from Bio.Blast import NCBIWWW

# configuration
Entrez.email = "curiousgeorge@ufl.edu"

# prompt user for organism and number of results desired
# query = input("Enter organism name: ")
# usermax = input("Enter maximum number of results: ")
# gene = input("Enter gene of interest: ")

# WHAT DO:
# User creates a file the specifies what sequences we're going to send/search in genbank
# Look into XML file format for genbank location tags

# TODO: For development CAN GET RID OF
query = "cantharellus minor"
usermax = 100
gene = "TEF"

# send query to Entrez nucleotide database
handle = Entrez.esearch(db="nucleotide", term=query, field="organism", retmax=usermax)
record = Entrez.read(handle)
handle.close()

# get the 'IDList' field from the results
ids = record["IdList"]

# function to get sequences from the genbank IDs returned from the search and write to 'query.fastas'
def get_sequences(ids):
    with open("query.fasta", "w") as out_file:
        for seq_id in ids:
            handle = Entrez.efetch(
                db="nucleotide", id=seq_id, rettype="fasta", retmode="text", retmax=1
            )
            out_file.write(handle.read())  # add data to end of file
    handle.close()


get_sequences(ids)

# function to get genbank data from the genbank IDs returned from the search and write to 'query_gb.txt'
def get_gb(ids):
    with open("query_gb.txt", "w") as out_file:
        for seq_id in ids:
            handle = Entrez.efetch(
                db="nucleotide", id=seq_id, rettype="genbank", retmode="text", retmax=1
            )
            out_file.write(handle.read())  # add data to end of file
    handle.close()


get_gb(ids)

# PART 2 - parsing the metadata from the sequences returned from our search, printing matching records...
# my_seqlist = []
# for seq_record in SeqIO.parse("query.fasta", "fasta"):
#     my_seqlist.append(seq_record)
my_seqlist = [seq_record for seq_record in SeqIO.parse("query.fasta", "fasta")]

# make new fasta file for seqs matching that gene
with open("second.fasta", "w") as out_file:
    for seq_record in my_seqlist:
        if gene in seq_record.description:
            out_file.write(">" + str(seq_record.id) + "\n")
            out_file.write(str(seq_record.seq) + "\n")


# Code adapted from Biopython help forums and adapted from work in the Biopython Cookbook Change et al 2020

# FIXME: NONE OF THIS WORKS?

# runs a loop to check if there are any duplicate sequences in the file
def remove_duplicates(records):
    check_sequences = set()
    for record in records:
        checksum = seguid(record.seq)
        if checksum in check_sequences:
            print(f"Removing duplicates {record.id} into outfile")
            continue
        check_sequences.add(checksum)
        yield record
        if record.id == record.seq:
            # if there are no duplicates a message will print no duplicates
            print("No duplicates!")


saved_sequences = remove_duplicates(SeqIO.parse("second.fasta", "fasta"))
count = SeqIO.write(saved_sequences, "final.fasta", "fasta")
print(f"Number of remaining sequences: {count}")

# read the new file and display the gene ids

# print("Unique designations are: ")

# final fasta will have non-duplicated sequences
for seq_record in SeqIO.parse("final.fasta", "fasta"):
    print(seq_record.id)

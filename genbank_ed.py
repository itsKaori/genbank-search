#!/usr/bin/env python
# coding: utf-8

# In[1]:


from Bio import Entrez, SeqIO
Entrez.email = 'curiousgeorge@ufl.edu' 

#prompt user for organism and number of results desired
query= input('Enter organism name ')
usermax=input('Enter maximum number of results ')

#send query to Entrez nucleotide database
handle = Entrez.esearch(db='nucleotide', term = query, field='organism', retmax=usermax)
record = Entrez.read(handle)
handle.close()

#get the 'IDList' field from the results
ids=(record["IdList"])

#function to get sequences from the genbank IDs returned from the search and write to 'query.fastas'
def get_sequences():
    open('query.fasta','w')
    for seq_id in ids:
        handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text",retmax=1)
        local_file=open('query.fasta','a') #a = append mode
        local_file.write(handle.read()) #add data to end of file
get_sequences()
handle.close()

#function to get genbank data from the genbank IDs returned from the search and write to 'query_gb.txt'
def get_gb():
    open('query.fasta','w')
    for seq_id in ids:
        handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="genbank", retmode="text",retmax=1)
        local_file=open('query_gb.txt','a') #a = append mode
        local_file.write(handle.read()) #add data to end of file
get_gb()
handle.close()


# In[ ]:



#In[ ]:

import sys #needed to be able to write a new file from standard output
#PART 2 - parsing the metadata from the sequences returned from our search, printing matching records... 
my_seqlist = []
for seq_record in SeqIO.parse("query.fasta", "fasta"):
    my_seqlist.append(seq_record)
    my_seqlist[0]
gene=input("Enter Gene of Interest: ")
sys.stdout = open("second.fasta", "w") #make new fasta file for seqs matching that gene
for seq_record in my_seqlist:
    if gene in seq_record.description:
        sys.stdout = open("second.fasta", "a")
        print ('>' +str(seq_record.id))
        print(seq_record.seq)
sys.stdout.close()

#In [ ]:

#Code adapted from Biopython help forums and adapted from work in the Biopython Cookbook Change et al 2020

#NONE OF THIS WORKS?

from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid

#runs a loop to check if there are any duplicate sequences in the file 
def remove_duplicates(records):
    check_sequences = set()
    for record in records:
        checksum = seguid(record.seq)
        if checksum in check_sequences:
            print("Removing duplicates %s into outfile" % record.id ) 
            continue
        check_sequences.add(checksum)
        yield record
        if record.id == record.seq:
            print("No duplicates!") #if there are no duplicates a message will print no duplicates

#saved_sequences = remove_duplicates(SeqIO.parse("second.fasta", "fasta")) 
#count = SeqIO.write(saved_sequences, "final.fasta", "fasta")
#print ("Number of remaining sequence %s" %count)

#read the new file and display the gene ids

#print("Unique designations are: ")

for seq_record in SeqIO.parse("final.fasta", "fasta"):   #final fasta will have non-duplicated sequences 
    print(seq_record.id)


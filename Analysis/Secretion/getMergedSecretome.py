# created by Gabriel on 26/01/16
# last modified by Gabriel on 26/01/16

# get one merged protein fasta file of secretome:

from Bio import SeqIO

path = "/path/to/Analysis/Secretion/"
with open(path + "MergedSecretome.fa", "w") as out_file:
    seq_all = []
    species = ("Mepe", "Usbr", "Usho", "Usma", "Spsc", "Sprs", "Sprz")
    for s in species:
        seqs = open(path + s + "/" + s + "_Secretome.fa")
        records = SeqIO.parse(seqs, "fasta")
        for record in records:
            seq_all.append(record)
    for i in range(len(seq_all)):
        out_file.write(">" + seq_all[i].id + "\n" + str(seq_all[i].seq) + "\n")

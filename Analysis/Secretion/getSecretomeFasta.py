# created by Gabriel on 22/02/17
# last modified by Gabriel on 22/02/17

# Get the fasta sequences for predicted secreted proteins

from Bio import SeqIO

species = ("Mepe", "Usho", "Usbr", "Usma", "Spsc", "Sprs", "Sprz")
for s in species:
    a = []
    with open("/path/to/Analysis/Secretion/" + s + "/" + s + "_Secretome.fa", "w") as out_file:
        f = open("/path/to/Analysis/Secretion/" + s + "/" + s + "_Secretome.txt")
        line = f.readlines()
        for l in line:
            a.append(l.rstrip())
        seqs = open("/path/to/Data/" + s + "/" + s + "_valid_prot.fa")
        records = SeqIO.parse(seqs, "fasta")
        for record in records:
            if record.id in a:
                #if record.id not in error:
                out_file.write(">" + record.id + "\n" + str(record.seq) + "\n")

# created by Gabriel on 22/02/17
# last modified by Gabriel on 22/02/17

# Get the fasta sequences for predicted secreted proteins

from Bio import SeqIO


#error = ("Mepe_4350", "Spre_5502") # I realized that this proteins contain X which is 
                                   # not accepted in some programms. Therefore, it is 
                                   # excluded from here on (Gabriel, 06.03.17) 
                                   # for the moment, this is ignored!!!
species = ("Mepe", "Usho", "Usbr", "Usma", "Spsc", "Sprs", "Sprz")
for s in species:
    a = []
    with open("/home/schweizerg/Science-Data/Genomes/FungalGenomes6/Analysis/Secretion/" + s + "/" + s + "_Secretome.fa", "w") as out_file:
        f = open("/home/schweizerg/Science-Data/Genomes/FungalGenomes6/Analysis/Secretion/" + s + "/" + s + "_Secretome.txt")
        line = f.readlines()
        for l in line:
            a.append(l.rstrip())
        seqs = open("/home/schweizerg/Science-Data/Genomes/FungalGenomes6/Data/" + s + "/" + s + "_valid_prot.fa")
        records = SeqIO.parse(seqs, "fasta")
        for record in records:
            if record.id in a:
                #if record.id not in error:
                out_file.write(">" + record.id + "\n" + str(record.seq) + "\n")
# created by Gabriel on 09/02/17
# last modified by Gabriel on 09/02/17
# last modified by Gabriel on 10/07/17 (include novel Mpennsylvanicum data and SrS)

from Bio import SeqIO

species = ("Mepe", "Usho", "Usbr", "Usma", "Spsc", "Sprz", "Sprs")
AA = ("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "*")

# First, list all protein sequences that do not start with "M" or have non-terminal stop codons.
count1 = 0
count2 = 0
with open("/path/to/data/Data/InvalidSequences.csv", "w") as out_file:
    out_file.write("Species,Gene\n")
    for s in species:
        seq = open("/home/schweizerg/Science-Data/Genomes/FungalGenomes6/Data/" + s + "/" + s + "_prot_raw.fa")
        records = SeqIO.parse(seq, "fasta")
        for record in records:
            if (str(record.seq)[0]).upper() != "M":
                out_file.write(s + "," + record.id + "\n")
                count1 += 1 # count1 = 122
            elif "*" in str(record.seq)[:-1]:
                out_file.write(s + "," + record.id + "\n")
                count2 += 1 #  count2 = 0


# Check if there are sequences that contain no canonical amino acids.

countNonAA = 0
nonAA = []
for s in species:
    seq = open("/path/to/data/Data/" + s + "/" + s + "_prot_raw.fa")
    records = SeqIO.parse(seq, "fasta")
    for record in records:
        r = str(record.seq)
        for i in range(len(r)):
            if r[i] not in AA and r[i] != "X":
                print s, "   ", record.id, "   ", r[i]
                nonAA.append(record.id)
                countNonAA += 1 #  countNonAA = 0
                
# Change record.ids and write clean sequences to files

count_all = 0
for s in ("Usbr", "Usho", "Usma", "Spsc", "Sprz"):
    count_record = 0
    with open("/path/to/data/Data/" + s + "/" + s + "ChangedGeneIds.csv", "w") as out_file:
        out_file.write("GeneIdOld,GeneIdOldShort,GeneIdNew\n")
        seq = open("/path/to/data/Data/" + s + "/" + s + "_prot_raw.fa")
        records = SeqIO.parse(seq, "fasta")
        for record in records:
            count_record += 1
            old = record.id
            new = s + "_" + ("0000" + str(count_record))[-4:]
            out_file.write(old + "," + old.split("@")[1]  +  "," +  new + "\n")
    count_record2 = 0
    with open("/path/to/data/Data/" + s + "/" + s + "_valid_prot.fa", "w") as out_file:
        seq = open("/path/to/data/Data/" + s + "/" + s + "_prot_raw.fa")
        records2 = SeqIO.parse(seq, "fasta")
        for record2 in records2:
            count_record2 += 1
            new2 = s + "_" + ("0000" + str(count_record2))[-4:]
            if str(record2.seq)[0] == "M":
                count_all += 1
                if str(record2.seq)[-1:] == "*":
                    out_file.write(">" + new2 + "\n" + str(record2.seq)[:-1] + "\n")
                else:
                    out_file.write(">" + new2 + "\n" + str(record2.seq) + "\n")



for s in ("Mepe", "Sprs"):
    count_record = 0
    with open("/path/to/data/Data/" + s + "/" + s + "ChangedGeneIds.csv", "w") as out_file:
        out_file.write("GeneIdOld,GeneIdOldShort,GeneIdNew\n")
        seq = open("/path/to/data/Data/" + s + "/" + s + "_prot_raw.fa")
        records = SeqIO.parse(seq, "fasta")
        for record in records:
            count_record += 1
            old = record.id
            new = s + "_" + ("0000" + str(count_record))[-4:]
            out_file.write(old + "," + old  +  "," +  new + "\n")
    count_record2 = 0
    with open("/path/to/data/Data/" + s + "/" + s + "_valid_prot.fa", "w") as out_file:
        seq = open("/path/to/data/Data/" + s + "/" + s + "_prot_raw.fa")
        records2 = SeqIO.parse(seq, "fasta")
        for record2 in records2:
            count_record2 += 1
            new2 = s + "_" + ("0000" + str(count_record2))[-4:]
            if str(record2.seq)[0] == "M":
                count_all += 1
                if str(record2.seq)[-1:] == "*":
                    out_file.write(">" + new2 + "\n" + str(record2.seq)[:-1] + "\n")
                else:
                    out_file.write(">" + new2 + "\n" + str(record2.seq) + "\n")

# In total, 47927 sequences have been written to files.
# For each species: 
# Mepe: 6868; Usho: 7085; Usma: 6754; Usbr: 7225;
# Spsc: 6685; Sprz: 6656; Sprs: 6654

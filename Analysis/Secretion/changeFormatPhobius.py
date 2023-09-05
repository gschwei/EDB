# created by Gabriel on 16/12/15
# last modified by Gabriel on 16/12/15

# Here, I change the phobius output format in such a way that R can read it.

species = ("Mepe", "Usbr", "Usho", "Usma", "Spsc", "Sprs", "Sprz")

for s in species:
    with open("/path/to/Analysis/Secretion/" + s + "/" + s + "_phobius_ChangedFormat.csv", "w") as out_file:
        p = open("/path/to/Analysis/Secretion/"  + s + "/" + s + "_phobius.txt")
        lines = p.readlines()
        for line in lines:
            if not line[0:7] == "SEQENCE": # The first line is the header and therefore ignored.
                l = (line.rstrip()).split(" ")
                if len(l) == 27:
                    out_file.write(l[0] + "," + l[23] + "," + l[25] + "," + l[26] + "\n")
                elif len(l) == 26:
                    out_file.write(l[0] + "," + l[22] + "," + l[24] + "," + l[25] + "\n")
                else:
                    print s

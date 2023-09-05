# created by Gabriel on 02/08/17
# last modified by Gabriel on 02/08/17

# change the Interpro prediction output so that each entry has the same number of columns.

with open("/path/to/Analysis/Interpro/InterproScanAllOut_mod.txt", "w") as out_file:
    f= open("/path/to/Analysis/Interpro/InterproScanAllOut.txt")
    for line in f.readlines():
        x = ""
        l = line.rstrip()
        x = l.split("\t")
        out_file.write(x[0] + "@" + x[1] + "@" + x[2] + "@" + x[3] + "@" + x[4] + "@" + x[5] + "@" + x[6] + "@" + x[7] + "@" + x[8] + "@" + x[9] + "@" + x[10] + "\n")

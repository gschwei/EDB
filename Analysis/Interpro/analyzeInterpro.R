# created by Gabriel on 07/03/17
# last modified by Gabriel on 07/03/17

intp1<-read.table("/home/schweizerg/Science-Data/Genomes/FungalGenomes5/Analysis/Interpro/InterproScanAll.txt", stringsAsFactors=FALSE, quote="", sep="\t")
intp2<-read.table("/home/schweizerg/Science-Data/Genomes/FungalGenomes5/Analysis/Interpro/InterproScanAll2.txt", stringsAsFactors=FALSE, quote="", sep="\t")
intp<-rbind(intp1, intp2)
intp<-intp[, c(1,4,5,6)]
colnames(intp)<-c("GeneIdNew", "Predictor", "ProfileId", "ProfileDesc")

fam_comp<-read.table("/home/schweizerg/Science-Data/Genomes/FungalGenomes5/Analysis/Homology/OrthoMCL/TestSettings/RefFamilies/Parse2/GetFamilies/FamilyComposition.csv", stringsAsFactors=FALSE, sep=" ", header=TRUE)
fam_comp<-merge(fam_comp, intp, by="GeneIdNew", all.x=TRUE)

core<-subset(fam_comp, CoreFamily == "yes")

unique_core<-core[!duplicated(core[6:8]), ]
write.table(unique_core, file="/home/schweizerg/Science-Data/Genomes/FungalGenomes5/Analysis/Interpro/UniqueCore.txt", sep="%%", quote=FALSE, row.names=FALSE)

# The data in unique_core were manually analyzed to check which predicted domain is really functional. Now the result of this
# analysis is opend and combined with the other information.

unique_core_fun<-read.table("/home/schweizerg/Science-Data/Genomes/FungalGenomes5/Analysis/Interpro/UniqueCore.csv", stringsAsFactors=FALSE, sep=";", header=TRUE)
core<-merge(core, unique_core_fun, by="ProfileId", all.x=TRUE)

# now get families with members without domain, with domain and mixed.

fam<-character()
domain<-character()
for (f in unique(core$Family)) {
  fam<-c(fam, f)
  x<-numeric()
  sub<-subset(core, Family == f)
  for (m in unique(sub$GeneIdNew)) {
    sub_mem<-subset(sub, GeneIdNew == m)
    if (sum(sub_mem$Function == "yes", na.rm=TRUE) > 0) {
      x<-c(x, 1)
    } else {
      x<-c(x, 0)
    }
  }
  if (sum(x) == length(unique(sub$GeneIdNew))) {
    domain<-c(domain, "dom")
  } else if (sum(x) == 0) {
    domain<-c(domain, "nod")
  } else {
    domain<-c(domain, "mix")
  }
}
fam_dom<-data.frame(Family=fam, Domain=domain, stringsAsFactors=FALSE)
core<-merge(core, fam_dom, by="Family", all.x=TRUE)
write.table(core, file="/home/schweizerg/Science-Data/Genomes/FungalGenomes5/Analysis/Interpro/InterproInfoCore.txt", quote=FALSE, sep="|", row.names=FALSE)


# we have many families of mixed type. Therefore, I give a closer look at those. If half of the members agree on one domain type, we consider this
# family to be functional.
# Which families are those [analysis of this question was done manually. Read the result here:]

mix_fam_res<-read.table("/home/schweizerg/Science-Data/Genomes/FungalGenomes5/Analysis/Interpro/MixedFamilyResult.csv", stringsAsFactors=FALSE, sep=";", header=TRUE)

fam2<-character()
dom2<-character()
for (f in unique(core$Family)) {
  fam2<-c(fam2, f)
  if (subset(core, Family == f)$Domain[1] == "dom") {
    dom2<-c(dom2, "yes")
  } else if (subset(core, Family == f)$Domain[1] == "nod") {
    dom2<-c(dom2, "no")
  } else if (subset(core, Family == f)$Domain[1] == "mix") {
    if (subset(mix_fam_res, Family == f)$Domain == "dom") {
      dom2<-c(dom2, "yes")
    } else if (subset(mix_fam_res, Family == f)$Domain == "nod") {
      dom2<-c(dom2, "no")
    } else {
      dom2<-c(dom2, NA)
    }
  } else {
    dom2<-c(dom2, NA)
  }
}
df2<-data.frame(Family=fam2, Domain2=dom2)
core<-merge(core, df2, by="Family", all.x=TRUE)
write.csv(core, file="/home/schweizerg/Science-Data/Genomes/FungalGenomes5/Analysis/Interpro/CoreInfo.csv")
# some numbers:
length(unique(subset(core, Domain2 == "yes")$Family)) # 79
length(unique(subset(core, Domain2 == "no")$Family)) # 50
# known effectors with domain:
unique(subset(core, Domain2 == "yes" & !(is.na(EffectorName)))$EffectorName)
# "pep4"   "exg1-1" "exg1-2" "egl1"   "s15-3"  "hum2"   "hum3"   "rrm67"  "gas1"   "prb1"   "cmu1"   "pit4"   "suc2"   "rsp2"   "mrpl19"
unique(subset(core, Domain2 == "no" & !(is.na(EffectorName)))$EffectorName)
# "eff1-8"  "eff1-9"  "eff1-4"  "eff1-7"  "eff1-11" "eff1-10" "eff1-2"  "eff1-3"  "eff1-1"  "rep1"    "stp2"    "apB73"   "rsp3"    "rsp1"   "stp3"    "pep1"    "pit2"    "see1"    "stp1" 

# Nr prots with domain in each species
sum(substr(unique(subset(core, Domain2 == "yes")$GeneIdNew), 1, 4) == "Mepe") # 134
sum(substr(unique(subset(core, Domain2 == "yes")$GeneIdNew), 1, 4) == "Usbr") # 133
sum(substr(unique(subset(core, Domain2 == "yes")$GeneIdNew), 1, 4) == "Usho") # 141
sum(substr(unique(subset(core, Domain2 == "yes")$GeneIdNew), 1, 4) == "Usma") # 149
sum(substr(unique(subset(core, Domain2 == "yes")$GeneIdNew), 1, 4) == "Spsc") # 145
sum(substr(unique(subset(core, Domain2 == "yes")$GeneIdNew), 1, 4) == "Spre") # 142

# Nr prots without domain in each species
sum(substr(unique(subset(core, Domain2 == "no")$GeneIdNew), 1, 4) == "Mepe") # 65
sum(substr(unique(subset(core, Domain2 == "no")$GeneIdNew), 1, 4) == "Usbr") # 76
sum(substr(unique(subset(core, Domain2 == "no")$GeneIdNew), 1, 4) == "Usho") # 67
sum(substr(unique(subset(core, Domain2 == "no")$GeneIdNew), 1, 4) == "Usma") # 75
sum(substr(unique(subset(core, Domain2 == "no")$GeneIdNew), 1, 4) == "Spsc") # 84
sum(substr(unique(subset(core, Domain2 == "no")$GeneIdNew), 1, 4) == "Spre") # 86

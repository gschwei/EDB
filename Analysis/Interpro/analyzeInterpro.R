# created by Gabriel on 02/08/17
# last modified by Gabriel on 02/08/17

# Here, I parse the output of Interpro scan predictions to see which protein has a functional domain.

ips<-read.table("/path/to/Analysis/Interpro/InterproScanAllOut_mod.txt", stringsAsFactors=FALSE, sep="@", quote="")
ips2<-ips[, c(1,4,5,7,8)]
colnames(ips2)<-c("GeneIdNew", "PredSoftware", "ProfileId", "DomainStart", "DomainEnd")

fam<-read.table("/path/to/Analysis/Homology/OrthoMCL/FamilyComposition.csv", stringsAsFactors=FALSE, sep=",", header=TRUE)

core_fam<-subset(fam, NrSpecies == 7)$Family
gras_fam<-subset(fam, NrSpecies == 6 & Mepe == 0)$Family
nrow(subset(fam, NrSpecies == 2 & Sprz > 0 & Usma > 0)) # there are no "maize-specific" families.

library(seqinr)

mem<-character()
unit<-character()
for (i in 1:398) {
  x<-character()
  f<-character()
  y<-character()
  x<-x<-unlist(strsplit(paste("0000", i, sep=""), ""))
  y<-paste(x[(length(x)-2):length(x)], collapse="")
  f<-paste("FAM_", y, sep="")
  seqs<-read.fasta(paste("/path/to/Analysis/Homology/OrthoMCL/Family/", f, ".fa", sep=""))
  for (n in attr(seqs, "name")) {
    unit<-c(unit, f)
    mem<-c(mem, n)
  }
}
famcomp<-data.frame(Family=unit, GeneIdNew=mem, stringsAsFactors=FALSE)


id<-list()
for (s in c("Mepe", "Usbr", "Usho", "Usma", "Spsc", "Sprs", "Sprz")) {
  id[[s]]<-read.table(paste("/path/to/Data/", s, "/", s, "ChangedGeneIds.csv", sep=""), stringsAsFactors=FALSE, sep=",", header=TRUE)
}
id_all<-rbind(id$Mepe, id$Usbr, id$Usho, id$Usma, id$Spsc, id$Sprs, id$Sprz)
id_all<-id_all[, 2:3]
famcomp<-merge(famcomp, id_all, by="GeneIdNew", all.x=TRUE)

gs<-character()
co<-character()
u<-character()
for (k in unique(famcomp$Family)) {
  u<-c(u, k)
  if (k %in% gras_fam) {
    gs<-c(gs, "yes")
  } else {
    gs<-c(gs, "no")
  }
  if (k %in% core_fam) {
    co<-c(co, "yes")
  } else {
    co<-c(co, "no")
  }
}
df<-data.frame(Family=u, Core=co, GrasSpec=gs, stringsAsFactors=FALSE)
famcomp<-merge(famcomp, df, by="Family", all.x=TRUE)

famcomp<-merge(famcomp, ips2, by="GeneIdNew", all.x=TRUE)
ip1<-subset(famcomp, Core == "yes")$ProfileId
ip2<-subset(famcomp, GrasSpec == "yes")$ProfileId
ip_sel<-unique(c(ip1, ip2))

dom_info1<-read.table("/path/to/Analysis/Interpro/UniqueCore.csv", stringsAsFactors=FALSE, sep=";", header=TRUE)
ip_check<-ip_sel[!(ip_sel %in% dom_info1$ProfileId)]
write.table(ip_check, file="/path/to/Analysis/Interpro/DomainCheck.txt", quote=FALSE, row.names=FALSE)

dom_info2<-read.table("/path/to/Analysis/Interpro/DomainDecision.csv", stringsAsFactors=FALSE, sep=";", header=TRUE)
colnames(dom_info2)<-c("ProfileId", "Function")
dom_info<-rbind(dom_info1, dom_info2)

candidates<-subset(famcomp, Core == "yes" | GrasSpec == "yes")

func<-character()
for (j in 1:nrow(candidates)) {
  if (is.na(candidates[j, "ProfileId"])) {
    func<-c(func, "no")
  } else {
    func<-c(func, subset(dom_info, ProfileId == candidates[j, "ProfileId"])$Function)
  }
}
candidates["Function"]<-func

sq<-character()
dom<-character()

for (d in unique(candidates$GeneIdNew)) {
  sq<-c(sq, d)
  sub<-subset(candidates, GeneIdNew == d)
  if ("yes" %in% sub$Function) {
    dom<-c(dom, "yes")
  } else {
    dom<-c(dom, "no")
  }
}
df_dom<-data.frame(GeneIdNew=sq, Function2=dom, stringsAsFactors=FALSE)

candidates<-merge(candidates, df_dom, by="GeneIdNew", all.x=TRUE)

gr<-character()
dt<-character()
for (e in unique(candidates$Family)) {
  gr<-c(gr, e)
  sub<-subset(candidates, Family == e)
  if (sum(sub$Function2 == "yes") == nrow(sub)) {
    dt<-c(dt, "Domain")
  } else if (sum(sub$Function2 == "no") == nrow(sub)) {
    dt<-c(dt, "NoDomain")
  } else {
    dt<-c(dt, "Mix")
  }
}

df_dt<-data.frame(Family=gr, Domain=dt)

candidates<-merge(candidates, df_dt, by="Family", all.x=TRUE)


# get a table only with U. maydis entries:
um_cand<-subset(candidates, substr(GeneIdNew, 1, 4) == "Usma")

umf<-character()
umdt<-character()
for (e in unique(um_cand$Family)) {
  umf<-c(umf, e)
  sub<-subset(um_cand, Family == e)
  if (sum(sub$Function2 == "yes") == nrow(sub)) {
    umdt<-c(umdt, "Domain")
  } else if (sum(sub$Function2 == "no") == nrow(sub)) {
    umdt<-c(umdt, "NoDomain")
  } else {
    umdt<-c(umdt, "Mix")
  }
}
um_df_dt<-data.frame(Family=umf, Domain2=umdt, stringsAsFactors=FALSE)

um_cand<-merge(um_cand, um_df_dt, by="Family", all.x=TRUE)

um<-character()
cl<-character()
umn<-numeric()
umd<-character()
umt<-character()

for (p in unique(um_cand$Family)) {
  cl<-c(cl, p)
  um_sub<-subset(um_cand, Family == p)
  um<-c(um, paste(unique(um_sub$GeneIdOldShort), collapse=";"))
  umn<-c(umn, length(unique(um_sub$GeneIdOldShort)))
  umd<-c(umd, unique(um_sub$Domain2))
  if (unique(um_sub$Core) == "yes") {
    umt<-c(umt, "Core")
  }
  if (unique(um_sub$GrasSpec) == "yes") {
    umt<-c(umt, "GrassSpec")
  }
}
um_tbl<-data.frame(Family=cl, NrUmGenes=umn, UmGenes=um, DomainType=umd, Selection=umt, stringsAsFactors=FALSE)
write.table(um_tbl, file="/path/to/Analysis/Interpro/UmSummary.txt", quote=FALSE, row.names=FALSE, sep=",")

# Now I add to each gene ID old short in candidates the information whether the gene has a functional domain or not.

mt<-character()
for (k in 1:nrow(candidates)) {
  if (candidates[k, "Function2"] == "yes") {
    mt<-c(mt, paste(candidates[k, "GeneIdOldShort"], "@wd", sep=""))
  } else {
    mt<-c(mt, paste(candidates[k, "GeneIdOldShort"], "@nd", sep=""))
  }
}
candidates["GeneIdOldShortDom"]<-mt
write.table(candidates, file="/path/to/Analysis/Interpro/Candidates.txt", quote=FALSE, sep=";", row.names=FALSE)


xx<-unique(subset(candidates, Core == "yes")$Family)
yy<-unique(subset(candidates, GrasSpec == "yes")$Family)

sel<-character()
grop<-character()
for (i in 1:398) {
  x<-character()
  f<-character()
  y<-character()
  x<-x<-unlist(strsplit(paste("0000", i, sep=""), ""))
  y<-paste(x[(length(x)-2):length(x)], collapse="")
  f<-paste("FAM_", y, sep="")
  grop<-c(grop, f)
  if (f %in% xx) {
    sel<-c(sel, "core")
  } else if (f %in% yy) {
    sel<-c(sel, "grass-specific")
  } else {
    sel<-c(sel, "accessory")
  }
}
df_zzz<-data.frame(Family=grop, Type=sel, stringsAsFactors=FALSE)
fc<-read.table("/path/to/Analysis/Homology/OrthoMCL/FamilyComposition.csv", stringsAsFactors=FALSE, sep=",", header=TRUE)
fc<-merge(fc, df_zzz, by="Family")

write.table(fc, file="/path/to/Analysis/Homology/OrthoMCL/GroupSummary.txt", quote=FALSE, row.names=FALSE)

# created by Gabriel on 22/02/17
# last modified by Gabriel on 11/07/17 # add Sprs as additional species

# Analyze the output of SignalP, TMHMM, Phobius and PSScan.

species<-"Usma"#c("Mepe", "Usbr", "Usho", "Usma", "Spsc", "Sprs", "Sprz")
for (spec in species) {
  tts<-numeric()
  tps<-numeric()
  s<-read.table(paste("/home/schweizerg/Science-Data/Genomes/FungalGenomes6/Analysis/Secretion/", spec, "/", spec, "_signalp.txt", sep=""), stringsAsFactors=FALSE)
  s<-s[, c(1,5,10)] # column 5 is Y-score
  colnames(s)<-c("Gene", "CS", "Desc")
  
  t<-read.table(paste("/home/schweizerg/Science-Data/Genomes/FungalGenomes6/Analysis/Secretion/", spec, "/", spec, "_TMHMM.txt", sep=""), stringsAsFactors=FALSE)
  t<-t[, c(1,5,6)]
  colnames(t)<-c("Gene", "PredHelT", "TopologyT")
  t["NrTMT"]<-as.integer(sapply(strsplit(t$PredHelT, "="), function(x) x[[2]]))
  e<-sapply(strsplit(t$TopologyT, "="), function(x) x[[2]])
  for (i in 1:nrow(t)) {
    x<-character()
    y<-character
    if (t[i, "NrTMT"] > 0) {
      a<-sapply(strsplit(e[i], "-"), function(k) k[[1]])
      if (length(unlist(a)) == 2) {
        tts<-c(tts, as.numeric(substr(a, 2, nchar(a))))
      } else {
        x<-unlist(strsplit(e[i], "-"))
        if ("o" %in% unlist(strsplit(x[length(x) -1], ""))) {
          y<-strsplit(x[length(x)-1], "o")
          tts<-c(tts, sapply(y, function(z) z[[2]]))
        } else {
          y<-strsplit(x[length(x)-1], "i")
          tts<-c(tts, sapply(y, function(z) z[[2]]))
        }
      }
    } else {
      tts<-c(tts, 0)
    }
  }
  t["TMTMaxStart"]<-as.integer(tts)
  t<-t[, c(1,5)]
  
  p<-read.table(paste("/home/schweizerg/Science-Data/Genomes/FungalGenomes6/Analysis/Secretion/", spec, "/", spec, "_phobius_ChangedFormat.csv", sep=""), stringsAsFactors=FALSE, sep=",")
  for (j in 1:nrow(p)) {
    if (p[j, "V2"] > 0 & p[j, "V3"] == 0) {
      m<-unlist(strsplit(p[j, "V4"], "-"))
      if (length(m) == 2) {
        b<-sapply(strsplit(p[j, "V4"], "-"), function(k) k[[1]])
        tps<-c(tps, as.numeric(substr(b, 2, nchar(b))))
      } else {
        n<-m[length(m) - 1]
        if ("o" %in% unlist(strsplit(n, ""))) {
          tps<-c(tps, as.numeric(sapply(strsplit(n, "o"), function(q) q[[2]])))
        } else {
          tps<-c(tps, as.numeric(sapply(strsplit(n, "i"), function(r) r[[2]])))
        }  
      }
    } else if (p[j, "V2"] > 0 & p[j, "V3"] != 0) {
      v<-unlist(strsplit(unlist(strsplit(p[j, "V4"], "/"))[2], "-"))
      w<-v[length(v)-1]
      if ("o" %in% unlist(strsplit(w, ""))) {
        tps<-c(tps, as.numeric(sapply(strsplit(w, "o"), function(g) g[[2]])))
      } else {
        tps<-c(tps, as.numeric(sapply(strsplit(w, "i"), function(h) h[[2]])))
      }
    } else {
      tps<-c(tps, 0)
    }
  }
  colnames(p)<-c("Gene", "TMP", "SPP", "DomainP")
  p["TMPMaxStart"]<-as.integer(tps)
  p<-p[, c(1,5)]
  
  er_seqs<-read.table(paste("/home/schweizerg/Science-Data/Genomes/FungalGenomes6/Analysis/Secretion/", spec, "/", spec, "_PSScan.txt", sep=""), stringsAsFactors=FALSE, sep=",")$V1
  
  summary<-merge(s, t, by="Gene")
  all<-merge(summary, p, by="Gene")
  sec<-subset(all, Desc == "Y" & TMTMaxStart < CS & TMPMaxStart < CS & !(Gene %in% er_seqs))
  rownames(sec)<-NULL
  write(sec$Gene, file=paste("/home/schweizerg/Science-Data/Genomes/FungalGenomes6/Analysis/Secretion/", spec, "/", spec, "_Secretome.txt", sep=""))
}
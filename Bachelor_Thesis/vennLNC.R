setwd('/home/aggelos/Desktop/Rotation/Files')

lnc <- read.delim("lncrna3", sep="\t", header=F)
names(lnc) <- c("gene", "id", "type")
head(lnc)


                                        #LNC P=0,001


NAFLvsNASHF2 <- read.delim("lnc_naflf2DEG0,001", sep="\t", header = T)
NAFLvsNASHF3 <- read.delim("lnc_naflf3DEG0,001", sep="\t", header = T)
NAFLvsNASHF4 <- read.delim("lnc_naflf4DEG0,001", sep="\t", header = T)

NAFLvsNASHF2  <- as.character(NAFLvsNASHF2$lnc_NAFLvsNASH.F2_DEGs)
NAFLvsNASHF3  <- as.character(NAFLvsNASHF3$lnc_NAFLvsNASH.F3_DEGs)
NAFLvsNASHF4  <- as.character(NAFLvsNASHF4$lnc_NAFLvsNASH.F4_DEGs)


F0F1vsF3 <- read.delim("lnc_f0f1vsf3DEG0,001", sep="\t", header = T)
F0F1vsF4 <- read.delim("lnc_f0f1vsf4DEG0,001", sep="\t", header = T)

F0F1vsF3  <- as.character(F0F1vsF3$lnc_NASH.F0.F1vsNASH.F3_DEGs)
F0F1vsF4  <- as.character(F0F1vsF4$lnc_NASH.F0.F1vsNASH.F4_DEGs)


common1 <-  Reduce(intersect, list(NAFLvsNASHF2,NAFLvsNASHF3,NAFLvsNASHF4))
common2 <-  Reduce(intersect, list(F0F1vsF3,F0F1vsF4))

common1 <- as.data.frame(common1)
common2 <- as.data.frame(common2)
names(common1) <- c("gene")
names(common2) <- c("gene")

mercommon1 <- merge(common1, lnc, by="gene")
mercommon2 <- merge(common2, lnc, by="gene")

write.table(mercommon1, "lnc_id_nafl0,001", sep="\t", col.names=T, row.names=T, quote=F)

write.table(mercommon2, "lnc_id_f0f1vsf340,001", sep="\t", col.names=T, row.names=T, quote=F)

                                        #LNC P=0,005


NAFLvsNASHF2 <- read.delim("lnc_naflf2DEG0,005", sep="\t", header = T)
NAFLvsNASHF3 <- read.delim("lnc_naflf3DEG0,005", sep="\t", header = T)
NAFLvsNASHF4 <- read.delim("lnc_naflf4DEG0,005", sep="\t", header = T)

NAFLvsNASHF2  <- as.character(NAFLvsNASHF2$lnc_NAFLvsNASH.F2_DEGs)
NAFLvsNASHF3  <- as.character(NAFLvsNASHF3$lnc_NAFLvsNASH.F3_DEGs)
NAFLvsNASHF4  <- as.character(NAFLvsNASHF4$lnc_NAFLvsNASH.F4_DEGs)

length(NAFLvsNASHF4)

F0F1vsF3 <- read.delim("lnc_f0f1vsf3DEG0,005", sep="\t", header = T)
F0F1vsF4 <- read.delim("lnc_f0f1vsf4DEG0,005", sep="\t", header = T)

F0F1vsF3  <- as.character(F0F1vsF3$lnc_NASH.F0.F1vsNASH.F3_DEGs)
F0F1vsF4  <- as.character(F0F1vsF4$lnc_NASH.F0.F1vsNASH.F4_DEGs)

length(F0F1vsF4)

common1 <-  Reduce(intersect, list(NAFLvsNASHF2,NAFLvsNASHF3,NAFLvsNASHF4))
common2 <-  Reduce(intersect, list(F0F1vsF3,F0F1vsF4))

common1 <- as.data.frame(common1)
common2 <- as.data.frame(common2)

names(common1) <- c("gene")
names(common2) <- c("gene")

mercommon1 <- merge(common1, lnc, by="gene")
mercommon2 <- merge(common2, lnc, by="gene")

write.table(mercommon1, "lnc_id_nafl0,005", sep="\t", col.names=T, row.names=T, quote=F)

write.table(mercommon2, "lnc_id_f0f1vsf340,005", sep="\t", col.names=T, row.names=T, quote=F)

                                        #LNC P=0,01


NAFLvsNASHF2 <- read.delim("lnc_naflf2DEG0,01", sep="\t", header = T)
NAFLvsNASHF3 <- read.delim("lnc_naflf3DEG0,01", sep="\t", header = T)
NAFLvsNASHF4 <- read.delim("lnc_naflf4DEG0,01", sep="\t", header = T)

NAFLvsNASHF2  <- as.character(NAFLvsNASHF2$lnc_NAFLvsNASH.F2_DEGs)
NAFLvsNASHF3  <- as.character(NAFLvsNASHF3$lnc_NAFLvsNASH.F3_DEGs)
NAFLvsNASHF4  <- as.character(NAFLvsNASHF4$lnc_NAFLvsNASH.F4_DEGs)

F0F1vsF3 <- read.delim("lnc_f0f1vsf3DEG0,01", sep="\t", header = T)
F0F1vsF4 <- read.delim("lnc_f0f1vsf4DEG0,01", sep="\t", header = T)

F0F1vsF3  <- as.character(F0F1vsF3$lnc_NASH.F0.F1vsNASH.F3_DEGs)
F0F1vsF4  <- as.character(F0F1vsF4$lnc_NASH.F0.F1vsNASH.F4_DEGs)

common1 <-  Reduce(intersect, list(NAFLvsNASHF2,NAFLvsNASHF3,NAFLvsNASHF4))
common2 <-  Reduce(intersect, list(F0F1vsF3,F0F1vsF4))

common1 <- as.data.frame(common1)
common2 <- as.data.frame(common2)
names(common1) <- c("gene")
names(common2) <- c("gene")

mercommon1 <- merge(common1, lnc, by="gene")
mercommon2 <- merge(common2, lnc, by="gene")

write.table(mercommon1, "lnc_id_nafl0,01", sep="\t", col.names=T, row.names=T, quote=F)

write.table(mercommon2, "lnc_id_f0f1vsf340,01", sep="\t", col.names=T, row.names=T, quote=F)


                                        #Intersect Common

nafl001 <- read.delim("lnc_id_nafl0,001", sep="\t" , header=T)
nafl005 <- read.delim("lnc_id_nafl0,005", sep="\t" , header=T)
nafl01 <- read.delim("lnc_id_nafl0,01", sep="\t" , header=T)

f0f1001 <- read.delim("lnc_id_f0f1vsf340,001", sep="\t" , header=T)
f0f1005 <- read.delim("lnc_id_f0f1vsf340,005", sep="\t" , header=T)
f0f101 <- read.delim("lnc_id_f0f1vsf340,01", sep="\t" , header=T)

common1  <- merge(nafl001, f0f1001)
common2  <- merge(nafl005, f0f1005)
common3  <- merge(nafl01, f0f101)

write.table(common3, "lnc_id_intersect0,01", sep="\t", col.names=T, row.names=T, quote=F)


                                        #VENN


library("ggVennDiagram")
library(ggplot2)
library(dplyr)


x1 <- list(NAFLvsNASHF2, NAFLvsNASHF3, NAFLvsNASHF4)

ggVennDiagram(x1, category.names = c("NAFLvsNASHF2", "NAFLvsNASHF3", "NAFLvsNASHF4")) + ggplot2::scale_fill_gradient(low = "blue",high = "yellow")

x2 <- list(F0F1vsF3, F0F1vsF4)

ggVennDiagram(x2, category.names = c("NASHF0F1vsF3", "NASHF0F1vsNASHF4")) + ggplot2::scale_fill_gradient(low = "blue",high = "yellow")

dev.copy(pdf, "LNC 0,01 NAFLvsNASH.pdf")

dev.copy(pdf, "LNC 0,01 F0F1vsF34.pdf")

dev.off()


naflinter001 <- read.delim("lnc_id_nafl0,001", sep="\t", header = T)
naflinter005 <- read.delim("lnc_id_nafl0,005", sep="\t", header = T)
naflinter01 <- read.delim("lnc_id_nafl0,01", sep="\t", header = T)

f0f1inter001 <- read.delim("lnc_id_f0f1vsf340,001", sep="\t", header = T)
f0f1inter005 <- read.delim("lnc_id_f0f1vsf340,005", sep="\t", header = T)
f0f1inter01 <- read.delim("lnc_id_f0f1vsf340,01", sep="\t", header = T)

naflinter001

naflinter001 <- as.character(naflinter001$gene)
naflinter005 <-  as.character(naflinter005$gene)
naflinter01 <-  as.character(naflinter01$gene)

f0f1inter001 <-  as.character(f0f1inter001$gene)
f0f1inter005 <-  as.character(f0f1inter005$gene)
f0f1inter01 <-  as.character(f0f1inter01$gene)

x1 <- list(naflinter01, f0f1inter01)

ggVennDiagram(x1, category.names = c("NAFL Intersect", "F0F1 Intersect")) + ggplot2::scale_fill_gradient(low = "blue",high = "yellow")

dev.copy(pdf, "LNC 0,01 Intersect.pdf")

dev.copy(pdf, "LNC 0,01 F0F1vsF34.pdf")

dev.off()

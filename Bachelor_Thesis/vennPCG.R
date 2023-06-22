setwd('/home/aggelos/Desktop/Rotation/Files')

library("ggVennDiagram")
library(ggplot2)
library(dplyr)

                                        #PCG NAFL with p=0.001

NAFLvsNASHF2 <- read.delim("naflf2DEG0,001", sep="\t", header = T)
NAFLvsNASHF3 <- read.delim("naflf3DEG0,001", sep="\t", header = T)
NAFLvsNASHF4 <- read.delim("naflf4DEG0,001", sep="\t", header = T)



                                        #make dataframe readable for Venn

NAFLvsNASHF2  <- as.character(NAFLvsNASHF2$NAFLvsNASH.F2_DEGs)
NAFLvsNASHF3  <- as.character(NAFLvsNASHF3$naflf3results...gene.)
NAFLvsNASHF4  <- as.character(NAFLvsNASHF4$naflf4results...gene.)


x1 <- list(NAFLvsNASHF2, NAFLvsNASHF3, NAFLvsNASHF4)

ggVennDiagram(x1, category.names = c("NAFLvsNASHF2", "NAFLvsNASHF3", "NAFLvsNASHF4")) + ggplot2::scale_fill_gradient(low = "blue",high = "yellow")


                                        #PCG NASHF0F1 with p=0.001


F0F1vsF2 <- read.delim("f0f1vsf2DEG0,001", sep="\t", header = T)
F0F1vsF3 <- read.delim("f0f1vsf3DEG0,001", sep="\t", header = T)
F0F1vsF4 <- read.delim("f0f1vsf4DEG0,001", sep="\t", header = T)

F0F1vsF2  <- as.character(F0F1vsF2$NASH.F0.F1vsNASH.F2_DEGs)
F0F1vsF3  <- as.character(F0F1vsF3$NASH.F0.F1vsNASH.F3_DEGs)
F0F1vsF4  <- as.character(F0F1vsF4$f0f1vsf4results...gene)

x2 <- list(F0F1vsF2, F0F1vsF3, F0F1vsF4)

ggVennDiagram(x2, category.names = c("NASHF0F1VSF2", "NASHF0F1vsF3", "NASHF0F1vsNASHF4")) + ggplot2::scale_fill_gradient(low = "blue",high = "yellow")


                                        #PCG NAFL with p=0.005


NAFLvsNASHF2 <- read.delim("naflf2DEG0,005", sep="\t", header = T)
NAFLvsNASHF3 <- read.delim("naflf3DEG0,005", sep="\t", header = T)
NAFLvsNASHF4 <- read.delim("naflf4DEG0,005", sep="\t", header = T)

NAFLvsNASHF2  <- as.character(NAFLvsNASHF2$NAFLvsNASH.F2_DEGs)
NAFLvsNASHF3  <- as.character(NAFLvsNASHF3$NAFLvsNASH.F3_DEGs)
NAFLvsNASHF4  <- as.character(NAFLvsNASHF4$NAFLvsNASH.F4_DEGs)

x1 <- list(NAFLvsNASHF2, NAFLvsNASHF3, NAFLvsNASHF4)

ggVennDiagram(x1, category.names = c("NAFLvsNASHF2", "NAFLvsNASHF3", "NAFLvsNASHF4")) + ggplot2::scale_fill_gradient(low = "blue",high = "yellow")


                                        #PCG NASHF0F1 with p=0.005

F0F1vsF2 <- read.delim("f0f1vsf2DEG0,005", sep="\t", header = T)
F0F1vsF3 <- read.delim("f0f1vsf3DEG0,005", sep="\t", header = T)
F0F1vsF4 <- read.delim("f0f1vsf4DEG0,005", sep="\t", header = T)

F0F1vsF2  <- as.character(F0F1vsF2$NASH.F0.F1vsNASH.F2_DEGs)
F0F1vsF3  <- as.character(F0F1vsF3$NASH.F0.F1vsNASH.F3_DEGs)
F0F1vsF4  <- as.character(F0F1vsF4$NASH.F0.F1vsNASH.F4_DEGs)

x2 <- list(F0F1vsF2, F0F1vsF3, F0F1vsF4)

ggVennDiagram(x2, category.names = c("NASHF0F1VSF2", "NASHF0F1vsF3", "NASHF0F1vsNASHF4")) + ggplot2::scale_fill_gradient(low = "blue",high = "yellow")

                                        #PCG NAFL with p=0.01

NAFLvsNASHF2 <- read.delim("naflf2DEG0,01", sep="\t", header = T)
NAFLvsNASHF3 <- read.delim("naflf3DEG0,01", sep="\t", header = T)
NAFLvsNASHF4 <- read.delim("naflf4DEG0,01", sep="\t", header = T)

NAFLvsNASHF2  <- as.character(NAFLvsNASHF2$NAFLvsNASH.F2_DEGs)
NAFLvsNASHF3  <- as.character(NAFLvsNASHF3$NAFLvsNASH.F3_DEGs)
NAFLvsNASHF4  <- as.character(NAFLvsNASHF4$NAFLvsNASH.F4_DEGs)


x1 <- list(NAFLvsNASHF2, NAFLvsNASHF3, NAFLvsNASHF4)

ggVennDiagram(x1, category.names = c("NAFLvsNASHF2", "NAFLvsNASHF3", "NAFLvsNASHF4")) + ggplot2::scale_fill_gradient(low = "blue",high = "yellow")


                                        #PCG NASHF0F1 with p=0.01

F0F1vsF2 <- read.delim("f0f1vsf2DEG0,01", sep="\t", header = T)
F0F1vsF3 <- read.delim("f0f1vsf3DEG0,01", sep="\t", header = T)
F0F1vsF4 <- read.delim("f0f1vsf4DEG0,01", sep="\t", header = T)


F0F1vsF2  <- as.character(F0F1vsF2$NASH.F0.F1vsNASH.F2_DEGs)
F0F1vsF3  <- as.character(F0F1vsF3$NASH.F0.F1vsNASH.F3_DEGs)
F0F1vsF4  <- as.character(F0F1vsF4$NASH.F0.F1vsNASH.F4_DEGs)

x2 <- list(F0F1vsF2, F0F1vsF3, F0F1vsF4)

ggVennDiagram(x2, category.names = c("NASHF0F1VSF2", "NASHF0F1vsF3", "NASHF0F1vsNASHF4")) + ggplot2::scale_fill_gradient(low = "blue",high = "yellow")

#######STEP TO FIND COMMON GENES BEFORE MERGING

common2 <-  Reduce(intersect, list(F0F1vsF2,F0F1vsF3,F0F1vsF4))

write.table(common2, "f0f1vsf234common0,01", sep="\t", col.names=T, row.names=T, quote=F)
length(common2)

                                        #Both Baseline Venn

idnafl001 <- read.delim("id_nafl0,001", sep="\t", header=T)
idnafl005 <- read.delim("id_nafl0,005", sep="\t", header=T)
idnafl01 <- read.delim("id_nafl0,01", sep="\t", header=T)

idf0v234001 <- read.delim("id_f0f1vsf2340,001", sep="\t", header=T)
idf0v234005 <- read.delim("id_f0f1vsf2340,005", sep="\t", header=T)
idf0v23401 <- read.delim("id_f0f1vsf2340,01", sep="\t", header=T)

idf0v34001 <- read.delim("id_f0f1vsf340,001", sep="\t", header=T)
idf0v34005 <- read.delim("id_f0f1vsf340,005", sep="\t", header=T)
idf0v3401 <- read.delim("id_f0f1vsf340,01", sep="\t", header=T)


idnafl001 <- as.character(idnafl001$gene)
idnafl005 <- as.character(idnafl005$gene)
idnafl01 <- as.character(idnafl01$gene)

idf0v234001 <- as.character(idf0v234001$gene)
idf0v234005 <- as.character(idf0v234005$gene)
idf0v23401 <- as.character(idf0v23401$gene)

idf0v34001 <- as.character(idf0v34001$gene)
idf0v34005 <- as.character(idf0v34005$gene)
idf0v3401 <- as.character(idf0v3401$gene)

x1 <- list(idnafl01,idf0v3401)

ggVennDiagram(x1, category.names = c("NAFL Intersect", "F0F1 Intersect")) + ggplot2::scale_fill_gradient(low = "blue",high = "yellow")

dev.copy(pdf, "0,01 Intersect F34.pdf")

dev.off()

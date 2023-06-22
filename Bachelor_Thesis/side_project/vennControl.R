setwd('/home/aggelos/Desktop/Rotation/Files')

library("ggVennDiagram")
library(ggplot2)
library(dplyr)

ctrnafl <- read.delim("lnc_controlnaflresults0,01adj", sep="\t", header = T)
ctrf0f1 <- read.delim("lnc_controlf0f1results0,01adj", sep="\t", header = T)
ctrf2 <- read.delim("lnc_controlf2results0,01adj", sep="\t", header = T)
ctrf3 <- read.delim("lnc_controlf3results0,01adj", sep="\t", header = T)
ctrf4 <- read.delim("lnc_controlf4results0,01adj", sep="\t", header = T)


pctrnafl <- read.delim("controlnaflresults0,001adj", sep="\t", header = T)
pctrf0f1 <- read.delim("controlf0f1results0,001adj", sep="\t", header = T)
pctrf2 <- read.delim("controlf2results0,001adj", sep="\t", header = T)
pctrf3 <- read.delim("controlf3results0,001adj", sep="\t", header = T)
pctrf4 <- read.delim("controlf4results0,001adj", sep="\t", header = T)

#Find common DEGs from baseline Control (Intersect)

common1 <- as.data.frame(ctrnafl[,c("gene","id","type")], header=T)
common2 <- as.data.frame(ctrf0f1[,c("gene","id","type")], header=T)
common3 <- as.data.frame(ctrf2[,c("gene","id","type")], header=T)
common4 <- as.data.frame(ctrf3[,c("gene","id","type")], header=T)
common5 <- as.data.frame(ctrf4[,c("gene","id","type")], header=T)

common <- Reduce(intersect, list(common1, common2, common3, common4, common5))
nrow(common)

rownames(common) <- c(1:change) #depends on common

write.table(common, "id_control0,01", quote=F, sep="\t")

                                        #Make Dataframe for Venn

                                        #LNC

ctrlnafl  <- as.character(ctrnafl$gene)
ctrlf0f1  <- as.character(ctrf0f1$gene) 
ctrlf2  <- as.character(ctrf2$gene)
ctrlf3  <- as.character(ctrf3$gene)
ctrlf4  <- as.character(ctrf4$gene)

x1 <- list(ctrlnafl, ctrlf0f1, ctrlf2, ctrlf3, ctrlf4)

venn1 <- ggVennDiagram(x1, category.names = c("ControlvsNAFL", "ControlvsF0F1", "ControlvsF2","ControlvsF3","ControlvsF4")) + ggplot2::scale_fill_gradient(low = "blue",high = "yellow")
venn1

ggsave("LNC 0,01 Control.pdf", venn1)

                                        #PCG

pctrlnafl  <- as.character(pctrnafl$gene)
pctrlf0f1  <- as.character(pctrf0f1$gene) 
pctrlf2  <- as.character(pctrf2$gene)
pctrlf3  <- as.character(pctrf3$gene)
pctrlf4  <- as.character(pctrf4$gene)

x2 <- list(pctrlnafl, pctrlf0f1, pctrlf2, pctrlf3, pctrlf4)

venn2 <- ggVennDiagram(x2, category.names = c("ControlvsNAFL", "ControlvsF0F1", "ControlvsF2","ControlvsF3","ControlvsF4")) + ggplot2::scale_fill_gradient(low = "blue",high = "yellow")
venn2


ggsave("0,001 Control.pdf", venn2)


                                        #Extra Intersection testing

idcon01 <- read.delim("lnc_id_control0,01", sep = "\t", header=T) #177
idcon005 <- read.delim("lnc_id_control0,005", sep = "\t",header=T) #136 
idcon001 <- read.delim("lnc_id_control0,001", sep = "\t",header=T) #81
pidcon01 <- read.delim("id_control0,01", sep = "\t", header=T) #2625
pidcon005 <- read.delim("id_control0,005", sep = "\t", header=T) #2258
pidcon001 <- read.delim("id_control0,001", sep = "\t", header=T) #1566

                                        #Comparisons

idnafl01 <- read.delim("lnc_id_nafl0,01", sep = "\t", header=T) #9
idnafl005 <- read.delim("lnc_id_nafl0,005", sep = "\t", header=T) #6
idnafl001 <- read.delim("lnc_id_nafl0,001", sep = "\t", header=T) #1
pidnafl01 <- read.delim("id_nafl0,01", sep = "\t", header=T) #73
pidnafl005 <- read.delim("id_nafl0,005", sep = "\t", header=T)#46
pidnafl001 <-read.delim("id_nafl0,001", sep = "\t", header=T) #16

connafl01 <- merge(idcon01, idnafl01) #GAS6-AS1
connafl005 <- merge(idcon005, idnafl005) #GAS6-AS1
connafl001 <- merge(idcon001, idnafl001) #GAS6-AS1

pconnafl01 <- merge(pidcon01, pidnafl01) #ALDH3A1/ NAA / AKR1B10 / TTL / CCDC64 / C15orf52 / SULF2 / SLC22A12 / HPR
pconnafl005 <- merge(pidcon005, pidnafl005) #ALDH3A1 / NAA / AKR1B10
pconnafl001 <- merge(pidcon001, pidnafl001) #NAA


idf0f301 <- read.delim("lnc_id_f0f1vsf340,01", sep = "\t", header=T)#28
idf0f3005 <- read.delim("lnc_id_f0f1vsf340,005", sep = "\t", header=T) #17
idf0f3001 <- read.delim("lnc_id_f0f1vsf340,001", sep = "\t", header=T) #10

conf0f301 <- merge(idcon01, idf0f301)#GAS6-AS1
conf0f3005 <- merge(idcon005, idf0f3005)#GAS6-AS1
conf0f3001 <- merge(idcon001, idf0f3001)#GAS6-AS1


pidf0f201 <- read.delim("id_f0f1vsf2340,01", sep = "\t", header=T)
pidf0f2005 <- read.delim("id_f0f1vsf2340,005", sep = "\t", header=T)
pidf0f2001 <- read.delim("id_f0f1vsf2340,001", sep = "\t", header=T)

pidf0f301 <- read.delim("id_f0f1vsf340,01", sep = "\t", header=T)
pidf0f3005 <- read.delim("id_f0f1vsf340,005", sep = "\t", header=T)
pidf0f3001 <- read.delim("id_f0f1vsf340,001", sep = "\t", header=T)

pconf0f201 <- merge(pidcon01, pidf0f201) #EEF1A2
pconf0f2005 <- merge(pidcon005, pidf0f2005) #EEF1A2
pconf0f2001 <- merge(pidcon001, pidf0f2001) #0

pconf0f301 <- merge(pidcon01, pidf0f301) #35 (AKR1B10)
pconf0f3005 <- merge(pidcon005, pidf0f3005)#18 (AKR1B10)
pconf0f3001 <- merge(pidcon001, pidf0f3001)#8 (AKR1B10)


idintersect01 <- read.delim("lnc_id_intersect0,01", sep = "\t", header=T)
idintersect005 <- read.delim("lnc_id_intersect0,005", sep = "\t", header=T)
idintersect001 <- read.delim("lnc_id_intersect0,001", sep = "\t", header=T)

conintersect01 <- merge(idcon01, idintersect01) #GAS6-AS1
conintersect005 <- merge(idcon005, idintersect005)#GAS6-AS1
conintersect001 <- merge(idcon001, idintersect001)#GAS6-AS1


pidintersectf201 <- read.delim("id_intersectwithf20,01", sep = "\t", header=T)
pidintersectf2005 <- read.delim("id_intersectwithf20,005", sep = "\t", header=T)
pidintersectf2001 <- read.delim("id_intersectwithf20,001", sep = "\t", header=T)

pidintersect01 <- read.delim("id_intersectnof20,01", sep = "\t", header=T)
pidintersect005 <- read.delim("id_intersectnof20,005", sep = "\t", header=T)
pidintersect001 <- read.delim("id_intersectnof20,001", sep = "\t", header=T)

pconintersectf201 <- merge(pidcon01, pidintersectf201) #0
pconintersectf2005 <- merge(pidcon005, pidintersectf2005) #0
pconintersectf2001 <- merge(pidcon001, pidintersectf2001) #0

pconintersect01 <- merge(pidcon01, pidintersect01) #TTL / SULF2 / AKR1B10
pconintersect005 <- merge(pidcon005, pidintersect005) # AKR1B10
pconintersect001 <- merge(pidcon001, pidintersect001) # 0 



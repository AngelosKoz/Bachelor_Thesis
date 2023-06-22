setwd('/home/aggelos/Desktop/Rotation/Files')

library(dplyr)

conpat00 <- read.delim("conpat00", sep="\t", header = T) 
conpat <- read.delim("conpat", sep="\t", header = T)

prot <- read.delim("protein", sep="\t", header=T) #filtered pcgs (no 0 no NA)
lnc <- read.delim("lncrna", sep="\t", header = T) #filtered lnc (no 0 no NA)
lx <- read.delim("lncrna3", sep="\t", header=F)
px <- read.delim("protcod3", sep="\t", header = F)
names(lx) <- c("gene", "id", "type")
names(px) <- c("gene", "id", "type")

patstages <- read.delim("patstages", sep="\t", header=T)
patstagesorder <- read.delim("patorderstages", sep="\t", header=T)

prot <- prot[, -c(1,2,3)]
lnc <- lnc[, -c(1,2,3)]

head(prot)

mlprot <- as.data.frame(t(prot))

mllnc <- as.data.frame(t(lnc))

mlprot1 <- tibble::rownames_to_column(mlprot, "patient") #row names as column
mlprot2 <- mlprot1[ order(match(mlprot1$patient, patstagesorder$patient)), ] #

mllnc1 <- tibble::rownames_to_column(mllnc, "patient")
mllnc2 <- mllnc1[ order(match(mllnc1$patient, patstagesorder$patient)), ] 

View(mlprot2[1:216,1, drop=FALSE])

View(mllnc2[1:216,1, drop=FALSE])

rownames(mlprot2) <- mlprot2[,1]
rownames(mllnc2) <- mllnc2[,1]

mlprot2 <- mlprot2[,-1]
mllnc2 <- mllnc2[,-1]

order <- patstagesorder[,-1]

mlprot3 <- cbind(mlprot2, order)
mllnc3 <- cbind(mllnc2, order)

write.table(mlprot3, "MLproteindf", sep="\t", col.names=T, row.names=T, quote=F)

write.table(mllnc3, "MLlncdf", sep="\t", col.names=T, row.names=T, quote=F)


mlprotein <- read.delim("MLproteindf", sep="\t", header=T) ##MAIN DF FOR PCG

mllncrna <- read.delim("MLlncdf", sep="\t", header=T)##MAIN DF FOR LNC

View(mllncrna[1:216,10905, drop=FALSE])



##Remove low expression/var

mlvar <- read.delim("MLlncSTAGEnaflvf4", header=T, sep="\t")

stage <- as.data.frame(mlvar[,"group"], header=T)
names(stage) <- "group"

mlvar1 <- mlvar[, -10903] #remove stage column
mlvar2 <- as.data.frame(t(mlvar1)) 
test1 <- mlvar1[, colSums(mlvar1)!=0] #remove columns with 0
mlvar2$max <- apply(mlvar2,1,max) #add column with max
test2 <- mlvar2[mlvar2$max > 1 , ] #remove rows that have a max of 1
test2 <- test2[, -66] #remove the max column
test3 <- as.data.frame(t(test2)) #change df again
mlvarx <- cbind(test3, stage) #add stage info

colnames(mlvarx)
nrow(mlvarx)

write.table(mlvarx, "MLlncSTAGEnaflvf4noZERO", sep="\t", quote=F)

mlvar3 <- sapply(test3, function(x) sd(x, na.rm=T) / mean(x, na.rm=T) * 100) ##coeff var
mlvar3
length(mlvar3)

hist(mlvar3, breaks=7628/10)

boxplot(mlvar3)

dev.copy(pdf, "boxplot.pdf")

dev.off()



##alternative## mlvar2$var <- sapply(mlvar2, var) ##check again the var function

############


contgroup1 <- mllncrna

ncol(contgroup1)

stage <- as.data.frame(contgroup1[,"group"], header=T)
names(stage) <- "group"

stage

contgroup2 <- contgroup1[, -c(10903:10906)]



contgroup3 <- as.data.frame(t(contgroup2)) 
contgroup3$max <- apply(contgroup3,1,max) 
test2 <- contgroup3[contgroup3$max > 1 , ]
colnames(test2)

test2 <- test2[, -217] 
test3 <- as.data.frame(t(test2)) 
contgroupfin1 <- cbind(test3, stage) 

tt1 <- sapply(test3, function(x) sd(x, na.rm=T) / mean(x, na.rm=T) * 100) ##coeff var
length(tt1)

nrow(test2)

rownames(test2)

boxplot(tt1)

colnames(contgroupfin1)
rownames(contgroupfin1)


write.table(contgroupfin1, "MLlncMULTInozeroSTAGE", quote=F, sep= "\t")

##Remove specific gene

removal <- read.delim("MLpcgCTRvsF0F1cleanSTAGE", sep="\t", header=T)

ncol(removal)

nrow(removal)

removal1 <- removal[, -which(names(removal) %in% c("ENSG00000221869", "ENSG00000108176", "ENSG00000105894","ENSG00000066135"))] ##remove by ensembl id

colnames(removal1)

write.table(removal1, "MLpcgCTRvsF0F1cleanREMOVE", sep="\t", quote=F)



##RBP GUT paper

library("readxl")

rbp1 <- read_excel("RBPrawcounts.xlsx")
rbp2 <- read_excel("RBPDEG.xlsx")
rbp3 <- read_excel("RBPh_list.xlsx")
rbp4 <- read_excel("RBPKEGG.xlsx")

write.table(rbp3, "rbplist", sep="\t", quote=F, col.names=T, row.names=T)

rbpraw <- read.delim("rbpraw", sep="\t", header=T)
rbpdeg <- read.delim("rbpDEG", sep="\t", header=T)
rbpkegg <- read.delim("rbpKEGG", sep="\t", header=T)
rbplist <- read.delim("rbplist", sep="\t", header=T)

rbplist1 <- rbplist[, -c(3:26)]

names(rbplist1) <- c("gene", "id")

names(rbpraw)[1] <- "gene"

rbprawlist <- merge(rbplist1, rbpraw, by="gene")

write.table(rbprawlist, "rbprawlist", sep="\t", quote=F)

y1 <- read.delim("controlnaflresults0,01adj", sep="\t", header=T)
y2 <- read.delim("controlf0f1results0,01adj", sep="\t", header=T)
y3 <- read.delim("controlf2results0,01adj", sep="\t", header=T)
y4 <- read.delim("controlf3results0,01adj", sep="\t", header=T)
y5 <- read.delim("controlf4results0,01adj", sep="\t", header=T)

x1 <- read.delim("naflf2results0,01adj", sep="\t", header=T)
x2 <- read.delim("naflf3results0,01adj", sep="\t", header=T)
x3 <- read.delim("naflf4results0,01adj", sep="\t", header=T)

z1 <- read.delim("f0f1vsf2results0,01adj", sep="\t", header=T)
z2 <- read.delim("f0f1vsf3results0,01adj", sep="\t", header=T)
z3 <- read.delim("f0f1vsf4results0,01adj", sep="\t", header=T)

rbpy1 <- merge(rbplist1, y1, by="gene")
rbpy2 <- merge(rbplist1, y2, by="gene")
rbpy3 <- merge(rbplist1, y3, by="gene")
rbpy4 <- merge(rbplist1, y4, by="gene")
rbpy5 <- merge(rbplist1, y5, by="gene")

rbpx1 <- merge(rbplist1, x1, by="gene")
rbpx2 <- merge(rbplist1, x2, by="gene")
rbpx3 <- merge(rbplist1, x3, by="gene")

rbpz1 <- merge(rbplist1, z1, by="gene")
rbpz2 <- merge(rbplist1, z2, by="gene")
rbpz3 <- merge(rbplist1, z3, by="gene")

sign1 <- read.delim("id_nafl0,01", sep="\t", header=T)
sign2 <- read.delim("id_intersectnof20,01", sep="\t", header=T)
sign3 <- read.delim("id_f0f1vsf340,01", sep="\t", header=T)
sign4 <- read.delim("id_f0f1vsf2340,01", sep="\t", header=T)

signrbp1 <- merge(sign1, rbplist1, by="gene")
signrbp2 <- merge(sign2, rbplist1, by="gene")
signrbp3 <- merge(sign3, rbplist1, by="gene")
signrbp4 <- merge(sign4, rbplist1, by="gene")

nrow(sign2)

protwhole <- read.delim("protein", sep="\t", header=T)

protwhole1 <- read.delim("protcod3", sep="\t", header=F)


nrow(protwhole1)

names(protwhole) <- c("gene", "id", "type")

names(protwhole1) <- c("gene", "id", "type")


ppr <- merge(protwhole, rbplist1, by="gene")
ppr1 <- merge(protwhole1, rbplist1, by="gene")
ppr1 <- ppr1[, -4]
names(ppr1) <- c("gene", "id", "type")
prot2 <- tibble::rownames_to_column(prot, "gene")
rbpprot <- merge(prot2, rbplist1, by="gene")
rbpprot1 <- rbpprot[, -218]
rbpprot2 <- rbpprot1[, -1]
rownames(rbpprot1) <- rbpprot1[,1]
rbpprot3 <- as.data.frame(t(rbpprot2))
rbpprot4 <- cbind(rbpprot3, patstagesorder)


rbpprotgroup <- rbpprot4[, -c(2171:2173, 2175)]

####################### RBP ############################

rbpgroup1 <- rbpprotgroup[-c(1:10, 62:202),]

stage <- as.data.frame(rbpgroup1[,"group"], header=T)
names(stage) <- "group"
rbpgroup2 <- rbpgroup1[, -2171] 
rbpgroup3 <- as.data.frame(t(rbpgroup2)) 
rbpgroup3$max <- apply(rbpgroup3,1,max) 
test2 <- rbpgroup3[rbpgroup3$max > 1 , ]
colnames(test2)
test2 <- test2[, -66] 
test3 <- as.data.frame(t(test2)) 
rbpgroupfin1 <- cbind(test3, stage) 

tt1 <- sapply(test3, function(x) sd(x, na.rm=T) / mean(x, na.rm=T) * 100) ##coeff var
tt1
length(tt1)

nrow(test2)

hist(tt1, breaks=2126/10)

boxplot(tt1)



#######RBP CONTROL

ncol(rbpprotgroup)


contgroup1 <- rbpprotgroup[-c(62:216),]

contgroup1 <- rbpprotgroup[-c(11:61,96:216),]

nrow(contgroup1)

stage <- as.data.frame(contgroup1[,"group"], header=T)
names(stage) <- "group"
stage

contgroup2 <- contgroup1[, -2171]


contgroup3 <- as.data.frame(t(contgroup2)) 
contgroup3$max <- apply(contgroup3,1,max) 
test2 <- contgroup3[contgroup3$max > 1 , ]
colnames(test2)

test2 <- test2[, -45] 
test3 <- as.data.frame(t(test2)) 
contgroupfin1 <- cbind(test3, stage) 

tt1 <- sapply(test3, function(x) sd(x, na.rm=T) / mean(x, na.rm=T) * 100) ##coeff var
length(tt1)

nrow(test2)

rownames(test2)

boxplot(tt1)

colnames(contgroupfin1)
rownames(contgroupfin1)


write.table(contgroupfin1, "MLrbpCTRvsF0F1cleanSTAGE", quote=F, sep= "\t")




                                        #PCG

mlproteinf0f4 <- mlprotein[-c(1:61, 96:202),]

mlproteinf0f4group1 <- mlproteinf0f4[, -c(18402:18403, 18405)]

colnames(mlproteinf0f4group1)

stage <- as.data.frame(mlproteinf0f4group1[,"group"], header=T)
names(stage) <- "group"



ncol(mlproteinnaflf4group1)

mlproteinnaflf4group2 <- mlproteinnaflf4group1[, -18402]

mlproteinnaflf4group3 <- as.data.frame(t(mlproteinnaflf4group2))

mlproteinnaflf4group3$max <- apply(mlproteinnaflf4group3,1,max)

test2 <- mlproteinnaflf4group3[mlproteinnaflf4group3$max > 1 , ]

colnames(test2)

test2 <- test2[, -49]

test3 <- as.data.frame(t(test2))

protgroupfin2 <- cbind(test3, stage) 

colnames(protgroupfin2)

write.table(protgroupfin2, "MLprotSTAGEf0f1vsf4STAGEclean", quote=F, sep="\t")

tt1 <- sapply(test3, function(x) sd(x, na.rm=T) / mean(x, na.rm=T) * 100)

tt1
length(tt1)


hist(tt1, breaks=17170/10)

boxplot(tt1)

## PCG Control


contgroup <- mlprotein[-c(62:216),]

contgroup1 <- contgroup[, -c(18402:18403, 18405)]

colnames(contgroup1)

stage <- as.data.frame(contgroup1[,"group"], header=T)
names(stage) <- "group"
stage



contgroup2 <- contgroup1[, -18402]


contgroup3 <- as.data.frame(t(contgroup2)) 
contgroup3$max <- apply(contgroup3,1,max) 
test2 <- contgroup3[contgroup3$max > 1 , ]
colnames(test2)

test2 <- test2[, -45] 
test3 <- as.data.frame(t(test2)) 
contgroupfin1 <- cbind(test3, stage) 

tt1 <- sapply(test3, function(x) sd(x, na.rm=T) / mean(x, na.rm=T) * 100) ##coeff var
length(tt1)

nrow(test2)

rownames(test2)

boxplot(tt1)

colnames(contgroupfin1)
rownames(contgroupfin1)


write.table(contgroupfin1, "MLrbpCTRvsNAFLcleanREMOVE", quote=F, sep= "\t")



setwd('/home/aggelos/Desktop/Rotation/Files')


pcg <- read.delim("protcod3", sep="\t", header=F)
lnc <- read.delim("lncrna3", sep="\t", header=F)

rbp <- read.delim("rbplist", sep="\t", header=T)

rbp <- rbp[, -c(3:26)]
names(rbp) <- c("gene", "id")


nrow(rbp)

rp11 <- pcg[grep("ENSG00000174243", pcg$gene),]
rp11

names(pcg) <- c("gene", "id", "type")
names(lnc) <- c("gene", "id", "type")



setwd('/home/aggelos/Desktop/Rotation/AML/Signatures')


x1 <- read.csv(file= "PCG_NAFL_PERF_signatures.csv", header=F)
x1 <- as.data.frame(x1[,-1])
names(x1) <- c("gene")
y1 <- merge(x1, pcg, by="gene")
y1

x2 <- read.csv(file= "PCG_F0F1vF4_PERF_REMOVE_signatures.csv", header=F)
x2 <- as.data.frame(x2[,-1])
names(x2) <- c("gene")
y2 <- merge(x2, pcg, by="gene")
y2

x21 <- read.csv(file= "PCG_F0F1_PERF_signatures.csv", header=F)
x21 <- as.data.frame(x21[,-1])
names(x21) <- c("gene")
y21 <- merge(x21, pcg, by="gene")
y21

x3 <- read.csv(file= "PCG_CONvNAFL_PERFsignatures.csv", header=F)
x3 <- as.data.frame(x3[,-1])
names(x3) <- c("gene")
y3 <- merge(x3, pcg, by="gene")
y3

x4 <- read.csv(file= "PCG_CONvF0_REMOVE_PERF_signatures.csv", header=F)
x4 <- as.data.frame(x4[,-1])
names(x4) <- c("gene")
y4 <- merge(x4, pcg, by="gene")
y4

x41 <- read.csv(file= "PCG_CONvF0_PERF_signatures.csv", header=F)
x41 <- as.data.frame(x41[,-1])
names(x41) <- c("gene")
y41 <- merge(x41, pcg, by="gene")
y41


x <- read.csv(file= ".csv", header=F)
x <- as.data.frame(x[,-1])
names(x) <- c("gene")
y <- merge(x, pcg, by="gene")
y

z1 <- merge(y2, y4)
z1


h1 <- read.delim("MLbothCTRLvsF0F1clean", sep="\t", header=T)

ncol(h1)

removal1 <- h1[, -which(names(h1) %in% c("ENSG00000066135","ENSG00000105894","ENSG00000108176","ENSG00000221869"))]

colnames(removal1)

write.table(removal1, "MLbothCTRLvsF0F1cleanREMOVE", sep="\t", quote=F)



#boxplot


c1 <- read.delim("MLpcgCTRvsNAFLcleanSTAGE", sep="\t", header=T)

rp11 <- pcg[grep("ENSG00000057657","ENSG00000074800","ENSG00000081041","ENSG00000128016","ENSG00000133639","ENSG00000141101","ENSG00000144043","ENSG00000153006","ENSG00000180389","ENSG00000197111","ENSG00000221869", pcg$gene),]



boxplot(c1[c1$group=="Control", grep("ENSG00000057657", colnames(c1))], c1[c1$group=="NAFL", grep("ENSG00000057657", colnames(c1))])


##################################################

f1 <- read.csv("BOTH_CONvF0_PERF_signatures.csv", sep=",", header = F)
f1 <- as.data.frame(f1[,-1])
names(f1) <- "gene"
f1

f1int <- read.csv("BOTH_CONvF0_INT_signatures.csv", sep=",", header = F)
f1int <- as.data.frame(f1int[,-1])
names(f1int) <- "gene"
f1int

f2 <- read.csv("BOTH_CONvF0_REMOVE_PERF_signatures.csv", sep=",", header = F)
f2 <- as.data.frame(f2[,-1])
names(f2) <- c("gene")
f2

f2int <- read.csv("BOTH_CONvF0_REMOVE_INT_signatures.csv", sep=",", header = F)
f2int <- as.data.frame(f2int[,-1])
names(f2int) <- c("gene")
f2int

f3 <- read.csv("BOTH_CONvNAFL_PERF_signatures.csv", sep=",", header = F)
f3 <- as.data.frame(f3[,-1])
names(f3) <- c("gene")
f3

f4 <- read.csv("BOTH_F0_REMOVE_PERF_signatures.csv", sep=",", header = F)
f4 <- as.data.frame(f4[,-1])
names(f4) <- c("gene")
f4

f4int <- read.csv("BOTH_F0_REMOVE_INT_signatures.csv", sep=",", header = F)
f4int <- as.data.frame(f4int[,-1])
names(f4int) <- c("gene")
f4int


f5 <- read.csv("BOTH_NAFLvF4_PERF_signatures.csv", sep=",", header = F)
f5 <- as.data.frame(f5[,-1])
names(f5) <- c("gene")
f5

f5int <- read.csv("BOTH_NAFLvF4_INT_signatures.csv", sep=",", header = F)
f5int <- as.data.frame(f5int[,-1])
names(f5int) <- c("gene")
f5int

g1 <- merge(f1 , pcg, by="gene") #con v f0
g1

g1int <- merge(f1int , pcg, by="gene") #con v f0
g1int


g2 <- merge(f2 , pcg, by="gene") #con v f0 removed
g2

g2l <- merge(f2, lnc, by="gene")
g2l


g3 <- merge(f3 , pcg, by="gene")#con v nafl
g3

g3l <- merge(f3, lnc, by="gene")
g3l

g4 <- merge(f4 , pcg, by="gene") #f0vf4 removed
g4

g4l <- merge(f4, lnc, by="gene")
g4l

g5 <- merge(f5 , pcg, by="gene")
g5

g5l <- merge(f5 , lnc, by="gene") #naflvf4 
g5l

g5int <- merge(f5int , pcg, by="gene")
g5int

g5intl <- merge(f5int, lnc, by="gene")
g5intl

rbptest <- merge(rbp, f2, by="gene")
rbptest


test <- merge(f2, f2int, by="gene")

##Overlaps pcg+lnc

o1 <- read.csv("PCG_CONvF0_PERF_signatures.csv", sep=",", header = F)
o1 <- as.data.frame(o1[,-1])
names(o1) <- "gene"
o1



o2 <- read.csv("PCG_CONvF0_REMOVE_PERF_signatures.csv", sep=",", header = F)
o2 <- as.data.frame(o2[,-1])
names(o2) <- "gene"
o2

o3 <- read.csv("PCG_CONvNAFL_PERFsignatures.csv", sep=",", header = F)
o3 <- as.data.frame(o3[,-1])
names(o3) <- "gene"
o3

o4 <- read.csv("PCG_F0F1_PERF_signatures.csv", sep=",", header = F)
o4 <- as.data.frame(o4[,-1])
names(o4) <- "gene"
o4


o5 <- read.csv("PCG_F0F1vF4_PERF_REMOVE_signatures.csv", sep=",", header = F)
o5 <- as.data.frame(o5[,-1])
names(o5) <- "gene"
o5

o6 <- read.csv("PCG_NAFL_PERF_signatures.csv", sep=",", header = F)
o6 <- as.data.frame(o6[,-1])
names(o6) <- "gene"
o6

o7 <- read.csv(".csv", sep=",", header = F)
o7 <- as.data.frame(o7[,-1])
names(o7) <- "gene"
o7

o8 <- read.csv(".csv", sep=",", header = F)
o8 <- as.data.frame(o8[,-1])
names(o8) <- "gene"
o8

o9 <- read.csv(".csv", sep=",", header = F)
o9 <- as.data.frame(o9[,-1])
names(o9) <- "gene"
o9

o10 <- read.csv(".csv", sep=",", header = F)
o10 <- as.data.frame(o10[,-1])
names(o10) <- "gene"
o10

o7 <- read.csv(".csv", sep=",", header = F)
o7 <- as.data.frame(o7[,-1])
names(o7) <- "gene"
o7

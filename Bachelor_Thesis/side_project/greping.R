setwd('/home/aggelos/Desktop/Rotation/Files')

y1 <- read.delim("controlnaflresultsoriginal", sep="\t", header=T)
y2 <- read.delim("controlf0f1resultsoriginal", sep="\t", header=T)
y3 <- read.delim("controlf2resultsoriginal", sep="\t", header=T)
y4 <- read.delim("controlf3resultsoriginal", sep="\t", header=T)
y5 <- read.delim("controlf4resultsoriginal", sep="\t", header=T)

ly1 <- read.delim("lnc_controlnaflresultsoriginal", sep="\t", header=T)
ly2 <- read.delim("lnc_controlf0f1resultsoriginal", sep="\t", header=T)
ly3 <- read.delim("lnc_controlf2resultsoriginal", sep="\t", header=T)
ly4 <- read.delim("lnc_controlf3resultsoriginal", sep="\t", header=T)
ly5 <- read.delim("lnc_controlf4resultsoriginal", sep="\t", header=T)

x1 <- read.delim("naflf0f1resultsoriginal", sep="\t", header=T)
x2 <- read.delim("naflf2resultsoriginal", sep="\t", header=T)
x3 <- read.delim("naflf3resultsoriginal", sep="\t", header=T)
x4 <- read.delim("naflf4resultsoriginal", sep="\t", header=T)

lx1 <- read.delim("lnc_naflf0f1resultsoriginal", sep="\t", header=T)
lx2 <- read.delim("lnc_naflf2resultsoriginal", sep="\t", header=T)
lx3 <- read.delim("lnc_naflf3resultsoriginal", sep="\t", header=T)
lx4 <- read.delim("lnc_naflf4resultsoriginal", sep="\t", header=T)

z1 <- read.delim("f0f1vsf2resultsoriginal", sep="\t", header=T)
z2 <- read.delim("f0f1vsf3resultsoriginal", sep="\t", header=T)
z3 <- read.delim("f0f1vsf4resultsoriginal", sep="\t", header=T)

lz1 <- read.delim("lnc_f0f1vsf2resultsoriginal", sep="\t", header=T)
lz2 <- read.delim("lnc_f0f1vsf3resultsoriginal", sep="\t", header=T)
lz3 <- read.delim("lnc_f0f1vsf4resultsoriginal", sep="\t", header=T)



cgas1 <- ly1[grep("GAS6-AS1", ly1$id),]
cgas2 <- ly2[grep("GAS6-AS1", ly2$id),]
cgas3 <- ly3[grep("GAS6-AS1", ly3$id),]
cgas4 <- ly4[grep("GAS6-AS1", ly4$id),]
cgas5 <- ly5[grep("GAS6-AS1", ly5$id),]

ngas1 <- lx1[grep("GAS6-AS1", lx1$id),]
ngas2 <- lx2[grep("GAS6-AS1", lx2$id),]
ngas3 <- lx3[grep("GAS6-AS1", lx3$id),]
ngas4 <- lx4[grep("GAS6-AS1", lx4$id),]


fgas1 <- lz1[grep("GAS6-AS1", lz1$id),]
fgas2 <- lz2[grep("GAS6-AS1", lz2$id),]
fgas3 <- lz3[grep("GAS6-AS1", lz3$id),]

head(cgas1) ##Find Gene ID 

cgas1 <- as.data.frame(cgas1[,c("log10padj","id","log2FoldChange")], header=T)
cgas2 <- as.data.frame(cgas2[,c("log10padj","id","log2FoldChange")], header=T)
cgas3 <- as.data.frame(cgas3[,c("log10padj","id","log2FoldChange")], header=T)
cgas4 <- as.data.frame(cgas4[,c("log10padj","id","log2FoldChange")], header=T)
cgas5 <- as.data.frame(cgas5[,c("log10padj","id","log2FoldChange")], header=T)


ngas1 <- as.data.frame(ngas1[,c("log10padj","id","log2FoldChange")], header=T)
ngas2 <- as.data.frame(ngas2[,c("log10padj","id","log2FoldChange")], header=T)
ngas3 <- as.data.frame(ngas3[,c("log10padj","id","log2FoldChange")], header=T)
ngas4 <- as.data.frame(ngas4[,c("log10padj","id","log2FoldChange")], header=T)


fgas1 <- as.data.frame(fgas1[,c("log10padj","id","log2FoldChange")], header=T)
fgas2 <- as.data.frame(fgas2[,c("log10padj","id","log2FoldChange")], header=T)
fgas3 <- as.data.frame(fgas3[,c("log10padj","id","log2FoldChange")], header=T)

head(cgas1) #Find id Column to change the name

names(cgas1)[2] <- "ENSG00000233695"
names(cgas2)[2] <- "ENSG00000233695"
names(cgas3)[2] <- "ENSG00000233695"
names(cgas4)[2] <- "ENSG00000233695"
names(cgas5)[2] <- "ENSG00000233695"

names(ngas1)[2] <- "ENSG00000233695"
names(ngas2)[2] <- "ENSG00000233695"
names(ngas3)[2] <- "ENSG00000233695"
names(ngas4)[2] <- "ENSG00000233695"

names(fgas1)[2] <- "ENSG00000233695"
names(fgas2)[2] <- "ENSG00000233695"
names(fgas3)[2] <- "ENSG00000233695"


row.names(cgas1) <- "CtrlvsNAFL"
row.names(cgas2) <- "CtrlvsF0F1"
row.names(cgas3) <- "CtrlvsF2"
row.names(cgas4) <- "CtrlvsF3"
row.names(cgas5) <- "CtrlvsF4"


row.names(ngas1) <- "NAFLvsF0F1"
row.names(ngas2) <- "NAFLvsF2"
row.names(ngas3) <- "NAFLvsF3"
row.names(ngas4) <- "NAFLvsF4"

row.names(fgas1) <- "F0F1vsF2"
row.names(fgas2) <- "F0F1vsF3"
row.names(fgas3) <- "F0F1vsF4"

allgas <- rbind(cgas1, cgas2, cgas3, cgas4, cgas5, ngas1, ngas2, ngas3, ngas4, fgas1, fgas2, fgas3)

allgas

write.table(allgas, "GAS-AS1_Info", quote=F, sep="\t")


##PAPER LINC

clinc1 <- ly1[grep("LINC00", ly1$id),]
clinc2 <- ly2[grep("LINC00", ly2$id),]
clinc3 <- ly3[grep("LINC00", ly3$id),]
clinc4 <- ly4[grep("LINC00", ly4$id),]
clinc5 <- ly5[grep("LINC00", ly5$id),]

nlinc1 <- lx1[grep("0326", lx1$id),]
nlinc2 <- lx2[grep("0326", lx2$id),]
nlinc3 <- lx3[grep("0326", lx3$id),]
nlinc4 <- lx4[grep("0326", lx4$id),]


flinc1 <- lz1[grep("LINC00326", lz1$id),]
flinc2 <- lz2[grep("LINC00326", lz2$id),]
flinc3 <- lz3[grep("LINC00326", lz3$id),]

head(nlinc1)

nrow(clinc1)

head(lz1)


ccct1 <- y1[grep("CCT3", y1$id),]
ccct2 <- y2[grep("CCT3", y2$id),]
ccct3 <- y3[grep("CCT3", y3$id),]
ccct4 <- y4[grep("CCT3", y4$id),]
ccct5 <- y5[grep("CCT3", y5$id),]

ncct1 <- x1[grep("CCT3", x1$id),]
ncct2 <- x2[grep("CCT3", x2$id),]
ncct3 <- x3[grep("CCT3", x3$id),]
ncct4 <- x4[grep("CCT3", x4$id),]


fcct1 <- z1[grep("CCT3", z1$id),]
fcct2 <- z2[grep("CCT3", z2$id),]
fcct3 <- z3[grep("CCT3", z3$id),]

head(ccct1)
head(ccct2)
head(ccct3)
head(ccct4)
head(ccct5)

head(ncct1)
head(ncct2)
head(ncct3)
head(ncct4)


head(fcct1)
head(fcct2)
head(fcct3)


##########  MLA ANALYSIS -- MIR  ################

cmir1 <- ly1[grep("MIR4435-1HG", ly1$id),]
cmir2 <- ly2[grep("MIR4435-1HG", ly2$id),]
cmir3 <- ly3[grep("MIR4435-1HG", ly3$id),]
cmir4 <- ly4[grep("MIR4435-1HG", ly4$id),]
cmir5 <- ly5[grep("MIR4435-1HG", ly5$id),]

nmir1 <- lx1[grep("MIR4435-1HG", lx1$id),]
nmir2 <- lx2[grep("MIR4435-1HG", lx2$id),]
nmir3 <- lx3[grep("MIR4435-1HG", lx3$id),]
nmir4 <- lx4[grep("MIR4435-1HG", lx4$id),]


fmir1 <- lz1[grep("MIR4435-1HG", lz1$id),]
fmir2 <- lz2[grep("MIR4435-1HG", lz2$id),]
fmir3 <- lz3[grep("MIR4435-1HG", lz3$id),]

head(cmir1) ##Find Gene ID 

cmir1 <- as.data.frame(cmir1[,c("log10padj","id","log2FoldChange")], header=T)
cmir2 <- as.data.frame(cmir2[,c("log10padj","id","log2FoldChange")], header=T)
cmir3 <- as.data.frame(cmir3[,c("log10padj","id","log2FoldChange")], header=T)
cmir4 <- as.data.frame(cmir4[,c("log10padj","id","log2FoldChange")], header=T)
cmir5 <- as.data.frame(cmir5[,c("log10padj","id","log2FoldChange")], header=T)


nmir1 <- as.data.frame(nmir1[,c("log10padj","id","log2FoldChange")], header=T)
nmir2 <- as.data.frame(nmir2[,c("log10padj","id","log2FoldChange")], header=T)
nmir3 <- as.data.frame(nmir3[,c("log10padj","id","log2FoldChange")], header=T)
nmir4 <- as.data.frame(nmir4[,c("log10padj","id","log2FoldChange")], header=T)


fmir1 <- as.data.frame(fmir1[,c("log10padj","id","log2FoldChange")], header=T)
fmir2 <- as.data.frame(fmir2[,c("log10padj","id","log2FoldChange")], header=T)
fmir3 <- as.data.frame(fmir3[,c("log10padj","id","log2FoldChange")], header=T)

head(cmir1) #Find id Column to change the name

names(cmir1)[2] <- "ENSG00000172965"
names(cmir2)[2] <- "ENSG00000172965"
names(cmir3)[2] <- "ENSG00000172965"
names(cmir4)[2] <- "ENSG00000172965"
names(cmir5)[2] <- "ENSG00000172965"

names(nmir1)[2] <- "ENSG00000172965"
names(nmir2)[2] <- "ENSG00000172965"
names(nmir3)[2] <- "ENSG00000172965"
names(nmir4)[2] <- "ENSG00000172965"

names(fmir1)[2] <- "ENSG00000172965"
names(fmir2)[2] <- "ENSG00000172965"
names(fmir3)[2] <- "ENSG00000172965"


row.names(cmir1) <- "CtrlvsNAFL"
row.names(cmir2) <- "CtrlvsF0F1"
row.names(cmir3) <- "CtrlvsF2"
row.names(cmir4) <- "CtrlvsF3"
row.names(cmir5) <- "CtrlvsF4"


row.names(nmir1) <- "NAFLvsF0F1"
row.names(nmir2) <- "NAFLvsF2"
row.names(nmir3) <- "NAFLvsF3"
row.names(nmir4) <- "NAFLvsF4"

row.names(fmir1) <- "F0F1vsF2"
row.names(fmir2) <- "F0F1vsF3"
row.names(fmir3) <- "F0F1vsF4"

allmir <- rbind(cmir1, cmir2, cmir3, cmir4, cmir5, nmir1, nmir2, nmir3, nmir4, fmir1, fmir2, fmir3)

allmir

write.table(allmir, "MIR4435-1HG_Info", quote=F, sep="\t")

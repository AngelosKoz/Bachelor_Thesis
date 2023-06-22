set.seed(1234)

setwd('/home/aggelos/Desktop/Rotation/Files')

pcg <- read.delim("protein", sep="\t", header=T)
lnc <- read.delim("lncrna", sep="\t", header=T)
rbp1 <- read.delim("rbplistID", sep="\t", header=T)

pint005naflr <- read.delim("controlnaflresults0,005adj", sep="\t", header=T)
pint001naflr <- read.delim("controlnaflresults0,001adj", sep="\t", header=T)
pint005f0f1r <- read.delim("controlf0f1results0,005adj", sep="\t", header=T)
pint001f0f1r <- read.delim("controlf0f1results0,001adj", sep="\t", header=T)

lint005naflr <- read.delim("lnc_controlnaflresults0,005adj", sep="\t", header=T)
lint001naflr <- read.delim("lnc_controlnaflresults0,001adj", sep="\t", header=T)
lint005f0f1r <- read.delim("lnc_controlf0f1results0,005adj", sep="\t", header=T)
lint001f0f1r <- read.delim("lnc_controlf0f1results0,001adj", sep="\t", header=T)

rint005naflr <- read.delim("RBP_controlnaflresults0,005adj", sep="\t", header=T)
rint001naflr <- read.delim("RBP_controlnaflresults0,001adj", sep="\t", header=T)
rint005f0f1r <- read.delim("RBP_controlf0f1results0,005adj", sep="\t", header=T)
rint001f0f1r <- read.delim("RBP_controlf0f1results0,001adj", sep="\t", header=T)


rbp <- merge(rbp1, pcg, by="id")
rbp <- rbp[,-2]
names(rbp)[2] <- "gene"
duplicated(rbp$gene)


pint005naflgene <- as.data.frame(pint005naflr[,1, drop=FALSE])
pint001naflgene <- as.data.frame(pint001naflr[,1, drop=FALSE])
pint005f0f1gene <- as.data.frame(pint005f0f1r[,1, drop=FALSE])
pint001f0f1gene <- as.data.frame(pint001f0f1r[,1, drop=FALSE])

lint005naflgene <- as.data.frame(lint005naflr[,1, drop=FALSE])
lint001naflgene <- as.data.frame(lint001naflr[,1, drop=FALSE])
lint005f0f1gene <- as.data.frame(lint005f0f1r[,1, drop=FALSE])
lint001f0f1gene <- as.data.frame(lint001f0f1r[,1, drop=FALSE])

rint005naflgene <- as.data.frame(rint005naflr[,1, drop=FALSE])
rint001naflgene <- as.data.frame(rint001naflr[,1, drop=FALSE])
rint005f0f1gene <- as.data.frame(rint005f0f1r[,1, drop=FALSE])
rint001f0f1gene <- as.data.frame(rint001f0f1r[,1, drop=FALSE])

pint005nafl <- merge(pcg, pint005naflgene, by="gene")
pint001nafl <- merge(pcg, pint001naflgene, by="gene")
pint005f0f1 <- merge(pcg, pint005f0f1gene, by="gene")
pint001f0f1 <- merge(pcg, pint001f0f1gene, by="gene")

lint005nafl <- merge(lnc, lint005naflgene, by="gene")
lint001nafl <- merge(lnc, lint001naflgene, by="gene")
lint005f0f1 <- merge(lnc, lint005f0f1gene, by="gene")
lint001f0f1 <- merge(lnc, lint001f0f1gene, by="gene")

rint005nafl <- merge(rbp, rint005naflgene, by="gene")
rint001nafl <- merge(rbp, rint001naflgene, by="gene")
rint005f0f1 <- merge(rbp, rint005f0f1gene, by="gene")
rint001f0f1 <- merge(rbp, rint001f0f1gene, by="gene")


pint005nafl1 <- pint005nafl[,c(4:219,1,2,3)] #re-order df
rownames(pint005nafl) <- pint005nafl[,1]

protlength <- read.delim("protlength2", sep="\t", header=F)
names(protlength) <- c("chr", "positionstart", "positionend", "strand", "gene")
rownames(protlength) <- protlength[,5]

protlength$length <- protlength$positionend - protlength$positionstart #add base count as column
print(protlength$chr)


library(DESeq2)
library(WGCNA)

options(stringAsFactors = FALSE) #Important setting!!
#data <- read.csv("LiverFemale3600.csv")

counts <- merge(pint005nafl1, protlength, by="gene", drop=FALSE)
rownames(counts) <- counts[,1]

rownames(metadata) <- metadata[,1]


dataExpr_deseq <- DESeqDataSetFromMatrix(countData=counts[,-c(1,218:224)], colData = metadata, design=~stage
mcols(dataExpr_deseq)$basepairs=counts$length

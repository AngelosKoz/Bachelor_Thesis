setwd('/home/aggelos/Desktop/Rotation/Files')

patients <- read.delim("conpat", sep="\t", header=T)

prot <- read.delim("protcod3", sep="\t", header=F)

names(prot) <- c("gene", "id", "type")

lnc <- read.delim("lncrna3", sep="\t", header=F)

names(lnc) <- c("gene", "id", "type")

patients <- patients[,-218]

protein <- merge(prot, patients, by="gene")

lncrna <- merge(lnc, patients,  by="gene")

row.names(protein) <- protein$gene

row.names(lncrna) <- lncrna$gene

protein <- protein[,-c(1,2,3)]

lncrna <- lncrna[,4:219]

library(DESeq2)

                                        #proteinDESeq2

proteinCounts <- protein

proteinframe <- data.frame(state=c(rep("control", 10), rep('patient', 206 )), replicate=c(1:10, 1:206), type=rep("single",216))

row.names(proteinframe) <-  names(proteinCounts)

proteinCol <- proteinframe[,c("state","type")]

proteindds <- DESeqDataSetFromMatrix(countData=proteinCounts, colData=proteinCol, design = ~ state)

colData(proteindds)$state <- factor(colData(proteindds)$state, levels=c("control", "patient"))

proteinDeseq <- DESeq(proteindds)

proteinresults <- results(proteinDeseq)

proteinresults <- data.frame(proteinresults)

proteinresults$gene <- row.names(proteinresults)

write.table(proteinresults, "proteinresults.ods",  sep="\t", col.names=T, row.names=T)

                                        #lncDESeq

lncCounts <- lncrna

lncframe <- data.frame(state=c(rep("control", 10), rep('patient', 206 )), replicate=c(1:10, 1:206), type=rep("single",216))

row.names(lncframe) <-  names(lncCounts)

lncCol <- lncframe[,c("state","type")]

lncdds <- DESeqDataSetFromMatrix(countData=lncCounts, colData=lncCol, design = ~ state)

colData(lncdds)$state <- factor(colData(lncdds)$state, levels=c("control", "patient"))

lncDeseq <- DESeq(lncdds)

lncresults <- results(lncDeseq)

lncresults <- data.frame(lncresults)

lncresults$gene <- row.names(lncresults)

write.table(lncresults, "lncresults.ods",  sep="\t", col.names=T, row.names=T)

                                        #removing 0 read counts


proteinresults <- read.delim("proteinresults.ods", sep="\t", header=T)

lncresults <- read.delim("lncresults.ods", sep="\t", header=T)

proteinresults0 <- proteinresults[!is.na(proteinresults$padj),]

lncresults0 <- lncresults[!is.na(lncresults$padj),]

proteinDEG <- as.data.frame(proteinresults0[,"gene"], header=T)

lncDEG <- as.data.frame(lncresults0[,"gene"], header=T)

names(proteinDEG)[1] <- "Pcg_Control_DEGs"

names(proteinDEG)[1] <- "Lnc_Control_DEGs"

patstages <- read.delim("patstages", sep="\t", header=T)

patstages$Patient <- gsub("C", "control", patstages$Patient)

patstages$Patient <- gsub("C", "control", patstages$Patient)

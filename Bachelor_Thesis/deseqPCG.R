setwd('/home/aggelos/Desktop/Rotation/Files')

patients <- read.delim("conpat", sep="\t", header=T)

prot <- read.delim("protcod3", sep="\t", header=F)
names(prot) <- c("gene", "id", "type")

protein <- merge(prot, patients, by="gene")

row.names(protein) <- protein$gene #add gene row name

protein <- protein[,-c(1,2,3)] #remove unwanted columns

nrow(protein)#check genes before remove = 21983

protein <- protein[ rowSums(protein)!=0, ] #remove 0 counts
nrow(protein) #check removed genes = 18401 // total of 3582 removed genes


                                        #divide based on disease stage


protein <- read.delim("protein", sep="\t", header=T) #read df with 0 counts
patstages <- read.delim("patstages", sep="\t", header=T)
controlstage <- subset(patstages, group == "Control", header=T) #extract specific stage /10
naflstage <- subset(patstages, group == "NAFL") # 51
nashf0f1stage <- subset(patstages, group == "NASH_F0-F1") # 34
nashf2stage <- subset(patstages, group == "NASH_F2") # 52
nashf3stage <- subset(patstages, group == "NASH_F3") # 55
nashf4stage <- subset(patstages, group == "NASH_F4") # 14

                                        #done manipulating second files/written tables

protein <- read.delim("protein", sep="\t", header=T)
patstages <- read.delim("patstages", sep="\t", header=T)
controlstage <- read.delim("controlstage", sep="\t", header=T)
naflstage <- read.delim("naflstage", sep="\t", header=T)
nashf0f1stage <- read.delim("nashf0f1stage", sep="\t", header=T)
nashf2stage <- read.delim("nashf2stage", sep="\t", header=T)
nashf3stage <- read.delim("nashf3stage", sep="\t", header=T)
nashf4stage <- read.delim("nashf4stage", sep="\t", header=T)

                                        #BASELINE NAFL


controlnafl <- rbind(controlstage, naflstage) #61
controlf0f1 <- rbind(controlstage, nashf0f1stage) #44
controlf2 <- rbind(controlstage, nashf2stage) #62
controlf3 <- rbind(controlstage, nashf3stage) #65
controlf4 <- rbind(controlstage, nashf4stage) #24

naflf0f1 <- rbind(naflstage, nashf0f1stage) #85
naflf2 <- rbind(naflstage, nashf2stage) #103
naflf3 <- rbind(naflstage, nashf3stage) #106
naflf4 <- rbind(naflstage, nashf4stage) #65

control_nafl <- protein[, as.character(names(protein)) %in% as.character(controlnafl$patient)]
control_f0f1 <- protein[, as.character(names(protein)) %in% as.character(controlf0f1$patient)]
control_f2 <- protein[, as.character(names(protein)) %in% as.character(controlf2$patient)]
control_f3 <- protein[, as.character(names(protein)) %in% as.character(controlf3$patient)]
control_f4 <- protein[, as.character(names(protein)) %in% as.character(controlf4$patient)]

nafl_f0f1 <- protein[, as.character(names(protein)) %in% as.character(naflf0f1$patient)] #removing the patients not included in our comparison from the whole gene count (protein)
nafl_f2 <- protein[, as.character(names(protein)) %in% as.character(naflf2$patient)]
nafl_f3 <- protein[, as.character(names(protein)) %in% as.character(naflf3$patient)]
nafl_f4 <- protein[, as.character(names(protein)) %in% as.character(naflf4$patient)]


                                        #Prepare variables for DESeq2
library(DESeq2)

controlnaflframe <- data.frame(stage=c(rep("control", 10), rep('nafl', 51)), replicate=c(1:10, 1:51), type=rep("single",61))
controlf0f1frame <- data.frame(stage=c(rep("control", 10), rep('nashf0f1', 34)), replicate=c(1:10, 1:34), type=rep("single",44))
controlf2frame <- data.frame(stage=c(rep("control", 10), rep('nashf2', 52)), replicate=c(1:10, 1:52), type=rep("single",62))
controlf3frame <- data.frame(stage=c(rep("control", 10), rep('nashf3', 55)), replicate=c(1:10, 1:55), type=rep("single",65))
controlf4frame <- data.frame(stage=c(rep("control", 10), rep('nashf4', 14)), replicate=c(1:10, 1:14), type=rep("single",24))

naflf0f1frame <- data.frame(stage=c(rep("nafl", 51), rep('nashf0f1', 34 )), replicate=c(1:51, 1:34), type=rep("single",85))
naflf2frame <- data.frame(stage=c(rep("nafl", 51), rep('nashf2', 52 )), replicate=c(1:51, 1:52), type=rep("single",103))
naflf3frame <- data.frame(stage=c(rep("nafl", 51), rep('nashf3', 55 )), replicate=c(1:51, 1:55), type=rep("single",106))
naflf4frame <- data.frame(stage=c(rep("nafl", 51), rep('nashf4', 14 )), replicate=c(1:51, 1:14), type=rep("single",65))

row.names(controlnaflframe) <- controlnafl[,1]
row.names(controlf0f1frame) <- controlf0f1[,1]
row.names(controlf2frame) <- controlf2[,1]
row.names(controlf3frame) <- controlf3[,1]
row.names(controlf4frame) <- controlf4[,1]

row.names(naflf0f1frame) <- naflf0f1[,1] #make the row names as the names of column 1
row.names(naflf2frame) <- naflf2[,1]
row.names(naflf3frame) <- naflf3[,1]
row.names(naflf4frame) <- naflf4[,1]


controlnaflCounts <- control_nafl[,match(rownames(controlnaflframe), colnames(control_nafl))]
controlf0f1Counts <- control_f0f1[,match(rownames(controlf0f1frame), colnames(control_f0f1))]
controlf2Counts <- control_f2[,match(rownames(controlf2frame), colnames(control_f2))]
controlf3Counts <- control_f3[,match(rownames(controlf3frame), colnames(control_f3))]
controlf4Counts <- control_f4[,match(rownames(controlf4frame), colnames(control_f4))]

naflf0f1Counts <- nafl_f0f1[,match(rownames(naflf0f1frame), colnames(nafl_f0f1))]#make the col names the same as row names
naflf2Counts <- nafl_f2[,match(rownames(naflf2frame), colnames(nafl_f2))]
naflf3Counts <- nafl_f3[,match(rownames(naflf3frame), colnames(nafl_f3))]
naflf4Counts <- nafl_f4[,match(rownames(naflf4frame), colnames(nafl_f4))]


                                        #deseq functions
controlnaflCol <- controlnaflframe[,c("stage","type")]
controlf0f1Col <- controlf0f1frame[,c("stage","type")]
controlf2Col <- controlf2frame[,c("stage","type")]
controlf3Col <- controlf3frame[,c("stage","type")]
controlf4Col <- controlf4frame[,c("stage","type")]

naflf0f1Col <- naflf0f1frame[,c("stage","type")]
naflf2Col <- naflf2frame[,c("stage","type")]
naflf3Col <- naflf3frame[,c("stage","type")]
naflf4Col <- naflf4frame[,c("stage","type")]

controlnafldds <- DESeqDataSetFromMatrix(countData=controlnaflCounts, colData=controlnaflCol, design = ~ stage)
controlf0f1dds <- DESeqDataSetFromMatrix(countData=controlf0f1Counts, colData=controlf0f1Col, design = ~ stage)
controlf2dds <- DESeqDataSetFromMatrix(countData=controlf2Counts, colData=controlf2Col, design = ~ stage)
controlf3dds <-DESeqDataSetFromMatrix(countData=controlf3Counts, colData=controlf3Col, design = ~ stage)
controlf4dds <- DESeqDataSetFromMatrix(countData=controlf4Counts, colData=controlf4Col, design = ~ stage)
    
naflf0f1dds <- DESeqDataSetFromMatrix(countData=naflf0f1Counts, colData=naflf0f1Col, design = ~ stage)
naflf2dds <- DESeqDataSetFromMatrix(countData=naflf2Counts, colData=naflf2Col, design = ~ stage)
naflf3dds <- DESeqDataSetFromMatrix(countData=naflf3Counts, colData=naflf3Col, design = ~ stage)
naflf4dds <- DESeqDataSetFromMatrix(countData=naflf4Counts, colData=naflf4Col, design = ~ stage)

colData(controlnafldds)$stage <- factor(colData(controlnafldds)$stage, levels=c("control", "nafl"))
colData(controlf0f1dds)$stage <- factor(colData(controlf0f1dds)$stage, levels=c("control", "nashf0f1"))
colData(controlf2dds)$stage <- factor(colData(controlf2dds)$stage, levels=c("control", "nashf2"))
colData(controlf3dds)$stage <- factor(colData(controlf3dds)$stage, levels=c("control", "nashf3"))
colData(controlf4dds)$stage <- factor(colData(controlf4dds)$stage, levels=c("control", "nashf4"))

colData(naflf0f1dds)$stage <- factor(colData(naflf0f1dds)$stage, levels=c("nafl", "nashf0f1"))
colData(naflf2dds)$stage <- factor(colData(naflf2dds)$stage, levels=c("nafl", "nashf2"))
colData(naflf3dds)$stage <- factor(colData(naflf3dds)$stage, levels=c("nafl", "nashf3"))
colData(naflf4dds)$stage <- factor(colData(naflf4dds)$stage, levels=c("nafl", "nashf4"))

###################Data Transformation for PCA###################

vsdT <- vst(controlf2dds, blind=TRUE) ##vst is faster

rld <- rlog(controlf2dds, blind=TRUE) ## if the size of the samples and therefore their size factors vary widely, the rlog transformation is a better option than VST

vsdF <- vst(controlf2dds, blind=FALSE)

rld <- rlog(controlf2dds, blind=FALSE)

head(assay(vsd), 3)

head(assay(rld), 3)

transvstFALSE <- assay(vsdF)

transvstTRUE <- assay(vsdT)

write.table(transvstFALSE, "transvstFALSE", sep="\t", quote=F)

write.table(transvstTRUE, "transvstTRUE", sep="\t", quote=F)

##END data for PCA

##Extra transformations



                                #Geometric Mean
library(psych)

gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){  ##using the = makes it a function
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}

controlnaflCounts$GeneGeoMeanControl <- apply(controlnaflCounts[,1:10],1,gm_mean)
controlnaflCounts$GeneGeoMeanNAFL <- apply(controlnaflCounts[,11:61],1, gm_mean)
controlf0f1Counts$GeneGeoMeanControl <- apply(controlf0f1Counts[,1:10],1,gm_mean)
controlf0f1Counts$GeneGeoMeanF0F1 <- apply(controlf0f1Counts[,11:44],1,gm_mean)
controlf2Counts$GeneGeoMeanControl <- apply(controlf2Counts[,1:10],1,gm_mean)
controlf2Counts$GeneGeoMeanF2 <- apply(controlf2Counts[,11:62],1,gm_mean)
controlf3Counts$GeneGeoMeanControl <- apply(controlf3Counts[,1:10],1,gm_mean)
controlf3Counts$GeneGeoMeanF3 <- apply(controlf3Counts[,11:65],1,gm_mean)
controlf4Counts$GeneGeoMeanControl <- apply(controlf4Counts[,1:10],1,gm_mean)
controlf4Counts$GeneGeoMeanF4 <- apply(controlf4Counts[,11:24],1,gm_mean)


ncol(controlnaflCounts)

x1 <- controlnaflCounts[,62:63]

x2 <- controlf0f1Counts[,45:46]

x3 <- controlf2Counts[,63:64]

x4 <- controlf3Counts[,66:67]

x5 <- controlf4Counts[,25:26]

xx <- cbind(x1,x2,x3,x4,x5)

xx <- xx[,-c(3,5,7,9)]

head(xx)

write.table(xx, "GeometricMeansPCG", sep="\t", quote=F)

y1 <- read.delim("controlnaflresultsoriginal", sep="\t", header=T)
y2 <- read.delim("controlf0f1resultsoriginal", sep="\t", header=T)
y3 <- read.delim("controlf2resultsoriginal", sep="\t", header=T)
y4 <- read.delim("controlf3resultsoriginal", sep="\t", header=T)
y5 <- read.delim("controlf4resultsoriginal", sep="\t", header=T)




controlnaflCounts$gene <- row.names(controlnaflCounts)
controlf0f1Counts$gene <- row.names(controlf0f1Counts)
controlf2Counts$gene <- row.names(controlf2Counts)
controlf3Counts$gene <- row.names(controlf3Counts)
controlf4Counts$gene <- row.names(controlf4Counts)

head(controlnaflCounts)

lncrna3 <- read.delim("lncrna", sep="\t", header=T)

lncrna3 <- lncrna3[, 1:3]


controlnaflresults <- merge(controlnaflresults, lncrna3)
controlf0f1results <- merge(controlf0f1results, lncrna3)
controlf2results <- merge(controlf2results, lncrna3)
controlf3results <- merge(controlf3results, lncrna3)
controlf4results <- merge(controlf4results, lncrna3)

controlnaflresults  <- merge(controlnaflresults, controlnaflCounts[,62:64], by="gene")
controlf0f1results  <- merge(controlf0f1results, controlf0f1Counts[,45:47], by="gene")
controlf2results  <- merge(controlf2results, controlf2Counts[,63:65], by="gene")
controlf3results  <- merge(controlf3results, controlf3Counts[,66:68], by="gene")
controlf4results  <- merge(controlf4results, controlf4Counts[,25:27], by="gene")

controlnaflresults$log10padj <- -log(controlnaflresults$padj,10)
controlf0f1results$log10padj <- -log(controlf0f1results$padj,10)
controlf2results$log10padj <- -log(controlf2results$padj,10)
controlf3results$log10padj <- -log(controlf3results$padj,10)
controlf4results$log10padj <- -log(controlf4results$padj,10)

row.names(controlnaflresults) <- controlnaflresults[,1]
row.names(controlf0f1results) <- controlf0f1results[,1]
row.names(controlf2results) <- controlf2results[,1]
row.names(controlf3results) <- controlf3results[,1]
row.names(controlf4results) <- controlf4results[,1]

##End extra

controlnaflDeseq <- DESeq(controlnafldds)

controlf0f1Deseq <- DESeq(controlf0f1dds)

controlf2Deseq <- DESeq(controlf2dds)

controlf3Deseq <- DESeq(controlf3dds)

controlf4Deseq <- DESeq(controlf4dds)

naflf0f1Deseq <- DESeq(naflf0f1dds)

naflf2Deseq <- DESeq(naflf2dds)

naflf3Deseq <- DESeq(naflf3dds)

naflf4Deseq <- DESeq(naflf4dds)

controlnaflresults <- results(controlnaflDeseq)
controlf0f1results <- results(controlf0f1Deseq)
controlf2results <- results(controlf2Deseq)
controlf3results <- results(controlf3Deseq)
controlf4results <- results(controlf4Deseq)

naflf0f1results <- results(naflf0f1Deseq)
naflf2results <- results(naflf2Deseq)
naflf3results <- results(naflf3Deseq)
naflf4results <- results(naflf4Deseq)

controlnaflresults <- data.frame(controlnaflresults)
controlf0f1results <- data.frame(controlf0f1results)
controlf2results <- data.frame(controlf2results)
controlf3results <- data.frame(controlf3results)
controlf4results <- data.frame(controlf4results)

naflf0f1results <- data.frame(naflf0f1results)
naflf2results <- data.frame(naflf2results)
naflf3results <- data.frame(naflf3results)
naflf4results <- data.frame(naflf4results)

controlnaflresults$gene <- row.names(controlnaflresults)
controlf0f1results$gene <- row.names(controlf0f1results)
controlf2results$gene <- row.names(controlf2results)
controlf3results$gene <- row.names(controlf3results)
controlf4results$gene <- row.names(controlf4results)

naflf0f1results$gene <- row.names(naflf0f1results) #add genes as row names,
naflf2results$gene <- row.names(naflf2results)
naflf3results$gene <- row.names(naflf3results)
naflf4results$gene <- row.names(naflf4results)


##Adding extra columns to original deseq

controlnaflCounts$GeneMeansControl <- apply(controlnaflCounts[,1:10],1,mean)
controlnaflCounts$GeneMeansNAFL <- apply(controlnaflCounts[,10:61],1,mean)
controlf0f1Counts$GeneMeansControl <- apply(controlf0f1Counts[,1:10],1,mean)
controlf0f1Counts$GeneMeansF0F1 <- apply(controlf0f1Counts[,10:44],1,mean)
controlf2Counts$GeneMeansControl <- apply(controlf2Counts[,1:10],1,mean)
controlf2Counts$GeneMeansF2 <- apply(controlf2Counts[,10:62],1,mean)
controlf3Counts$GeneMeansControl <- apply(controlf3Counts[,1:10],1,mean)
controlf3Counts$GeneMeansF3 <- apply(controlf3Counts[,10:65],1,mean)
controlf4Counts$GeneMeansControl <- apply(controlf4Counts[,1:10],1,mean)
controlf4Counts$GeneMeansF4 <- apply(controlf4Counts[,10:24],1,mean)


controlnaflCounts$gene <- row.names(controlnaflCounts)
controlf0f1Counts$gene <- row.names(controlf0f1Counts)
controlf2Counts$gene <- row.names(controlf2Counts)
controlf3Counts$gene <- row.names(controlf3Counts)
controlf4Counts$gene <- row.names(controlf4Counts)

head(controlnaflCounts)

head(prot)


controlnaflresults <- merge(controlnaflresults, prot)
controlf0f1results <- merge(controlf0f1results, prot)
controlf2results <- merge(controlf2results, prot)
controlf3results <- merge(controlf3results, prot)
controlf4results <- merge(controlf4results, prot)

controlnaflresults  <- merge(controlnaflresults, controlnaflCounts[,62:64], by="gene")
controlf0f1results  <- merge(controlf0f1results, controlf0f1Counts[,45:47], by="gene")
controlf2results  <- merge(controlf2results, controlf2Counts[,63:65], by="gene")
controlf3results  <- merge(controlf3results, controlf3Counts[,66:68], by="gene")
controlf4results  <- merge(controlf4results, controlf4Counts[,25:27], by="gene")

controlnaflresults$log10padj <- -log(controlnaflresults$padj,10)
controlf0f1results$log10padj <- -log(controlf0f1results$padj,10)
controlf2results$log10padj <- -log(controlf2results$padj,10)
controlf3results$log10padj <- -log(controlf3results$padj,10)
controlf4results$log10padj <- -log(controlf4results$padj,10)

row.names(controlnaflresults) <- controlnaflresults[,1]
row.names(controlf0f1results) <- controlf0f1results[,1]
row.names(controlf2results) <- controlf2results[,1]
row.names(controlf3results) <- controlf3results[,1]
row.names(controlf4results) <- controlf4results[,1]

write.table(controlnaflresults, "controlnaflresults0,001adj",  sep="\t", quote = F)

write.table(controlf0f1results, "controlf0f1results0,001adj",  sep="\t", quote = F)

write.table(controlf2results, "controlf2results0,001adj",  sep="\t", quote = F)

write.table(controlf3results, "controlf3results0,001adj",  sep="\t", quote = F)

write.table(controlf4results, "controlf4results0,001adj",  sep="\t", quote = F)

controlnaflresults <- read.delim("controlnaflresultsoriginal",  sep="\t", header=T)
controlf0f1results <- read.delim("controlf0f1resultsoriginal",  sep="\t", header=T)
controlf2results <- read.delim("controlf2resultsoriginal",  sep="\t", header=T)
controlf3results <- read.delim("controlf3resultsoriginal",  sep="\t", header=T)
controlf4results <- read.delim("controlf4resultsoriginal",  sep="\t", header=T)


nrow(controlf4results)

                                        #removing NA
controlnaflresults <- ctrnafl[!is.na(ctrnafl$padj),] #18401 initial, 16533 after
controlf0f1results <- ctrf0f1[!is.na(ctrf0f1$padj),] #18401 initial, 16404
controlf2results <- ctrf2[!is.na(ctrf2$padj),] #18401 initial, 16550
controlf3results <- ctrf3[!is.na(ctrf3$padj),] #18401 initial, 16614
controlf4results <- ctrf4[!is.na(ctrf4$padj),] #18401 initial, 15934

naflf0f1results <- naflf0f1results[!is.na(naflf0f1results$padj),] #18401 initial, 18014 after
naflf2results <- naflf2results[!is.na(naflf2results$padj),] #18401 initial, 16360 after
naflf3results <- naflf3results[!is.na(naflf3results$padj),] #18401 initial, 16413 after
naflf4results <- naflf4results[!is.na(naflf4results$padj),] #18401 initial, 16596 after

nrow(controlnaflresults)
nrow(controlf0f1results)
nrow(controlf2results)
nrow(controlf3results)
nrow(controlf4results)
                                 
                                        #Padj<0.01

controlnaflresults <- subset(controlnaflresults, padj<0.01) #5135
controlf0f1results <- subset(controlf0f1results, padj<0.01) #4771
controlf2results <- subset(controlf2results, padj<0.01) #5417
controlf3results <- subset(controlf3results, padj<0.01) #4854
controlf4results <- subset(controlf4results, padj<0.01) #4259


naflf0f1results <- subset(naflf0f1results, padj<0.01) #0 DEG vs 0
naflf2results <- subset(naflf2results, padj<0.01) #95 DEG vs 50
naflf3results <- subset(naflf3results, padj<0.01) #1791 DEG vs 907
naflf4results <- subset(naflf4results, padj<0.01) #1886 DEG vs 1369


                                        #Padj<0.005

controlnaflresults <- subset(controlnaflresults, padj<0.005) #4589
controlf0f1results <- subset(controlf0f1results, padj<0.005) #4238
controlf2results <- subset(controlf2results, padj<0.005) #4853
controlf3results <- subset(controlf3results, padj<0.005) #4275
controlf4results <- subset(controlf4results, padj<0.005) #3661

naflf0f1results <- subset(naflf0f1results, padj<0.005) #0 DEG vs 0
naflf2results <- subset(naflf2results, padj<0.005) #58 DEG vs 50
naflf3results <- subset(naflf3results, padj<0.005) #1389 DEG vs 907
naflf4results <- subset(naflf4results, padj<0.005) #1551 DEG vs 1369


                                        #removing values above padj=0.001

controlnaflresults <- subset(controlnaflresults, padj<0.001) #3672
controlf0f1results <- subset(controlf0f1results, padj<0.001) #3241
controlf2results <- subset(controlf2results, padj<0.001) #3815
controlf3results <- subset(controlf3results, padj<0.001) #3360
controlf4results <- subset(controlf4results, padj<0.001) #2555

naflf0f1results <- subset(naflf0f1results, padj<0.001) #0 DEG vs 0
naflf2results <- subset(naflf2results, padj<0.001) #21 DEG vs 50
naflf3results <- subset(naflf3results, padj<0.001) #915 DEG vs 907
naflf4results <- subset(naflf4results, padj<0.001) #1049 DEG vs 1369


                                        #finished results


naflf0f1results <- read.delim("naflf0f1resultsadj", sep="\t", header=T)
naflf2results <- read.delim("naflf2resultsadj", sep="\t", header=T)
naflf3results <- read.delim("naflf3resultsadj", sep="\t", header=T)
naflf4results <- read.delim("naflf4resultsadj", sep="\t", header=T)


                                        #extracting just the DEGs with 0.001 cutoff


naflf0f1DEG <- as.data.frame(naflf0f1results[,"gene"], header=T) #save specific column
naflf2DEG <- as.data.frame(naflf2results[,"gene"], header=T)
naflf3DEG <- as.data.frame(naflf3results[,"gene"], header=T)
naflf4DEG <- as.data.frame(naflf4results[,"gene"], header=T)

names(naflf0f1DEG)[1] <- "NAFLvsNASH-F0-F1_DEGs" ##rename column
names(naflf2DEG)[1] <- "NAFLvsNASH-F2_DEGs"
names(naflf3DEG)[1] <- "NAFLvsNASH-F3_DEGs"
names(naflf4DEG)[1] <- "NAFLvsNASH-F4_DEGs"



                                        #BASELINE NASHF0-F1 VS NASHF2-F4



f0f1vsf2 <- rbind(nashf0f1stage, nashf2stage) #86
f0f1vsf3 <- rbind(nashf0f1stage, nashf3stage) #89
f0f1vsf4 <- rbind(nashf0f1stage, nashf4stage) #48

f0f1_f2 <- protein[, as.character(names(protein)) %in% as.character(f0f1vsf2$patient)]
f0f1_f3 <- protein[, as.character(names(protein)) %in% as.character(f0f1vsf3$patient)]
f0f1_f4 <- protein[, as.character(names(protein)) %in% as.character(f0f1vsf4$patient)]

f0f1vsf2frame <- data.frame(stage=c(rep("nashf0f1", 34), rep('nashf2', 52 )), replicate=c(1:34, 1:52), type=rep("single",86))
f0f1vsf3frame <- data.frame(stage=c(rep("nashf0f1", 34), rep('nashf3', 55 )), replicate=c(1:34, 1:55), type=rep("single",89))
f0f1vsf4frame <- data.frame(stage=c(rep("nashf0f1", 34), rep('nashf4', 14 )), replicate=c(1:34, 1:14), type=rep("single",48))

row.names(f0f1vsf2frame) <- f0f1vsf2[,1]
row.names(f0f1vsf3frame) <- f0f1vsf3[,1]
row.names(f0f1vsf4frame) <- f0f1vsf4[,1]

f0f1vsf2Counts <- f0f1_f2[,match(rownames(f0f1vsf2frame), colnames(f0f1_f2))]
f0f1vsf3Counts <- f0f1_f3[,match(rownames(f0f1vsf3frame), colnames(f0f1_f3))]
f0f1vsf4Counts <- f0f1_f4[,match(rownames(f0f1vsf4frame), colnames(f0f1_f4))]

f0f1vsf2Col <- f0f1vsf2frame[,c("stage","type")]
f0f1vsf3Col <- f0f1vsf3frame[,c("stage","type")]
f0f1vsf4Col <- f0f1vsf4frame[,c("stage","type")]


f0f1vsf2dds <- DESeqDataSetFromMatrix(countData=f0f1vsf2Counts, colData=f0f1vsf2Col, design = ~ stage)
f0f1vsf3dds <- DESeqDataSetFromMatrix(countData=f0f1vsf3Counts, colData=f0f1vsf3Col, design = ~ stage)
f0f1vsf4dds <- DESeqDataSetFromMatrix(countData=f0f1vsf4Counts, colData=f0f1vsf4Col, design = ~ stage)



colData(f0f1vsf2dds)$stage <- factor(colData(f0f1vsf2dds)$stage, levels=c("nashf0f1", "nashf2"))
colData(f0f1vsf3dds)$stage <- factor(colData(f0f1vsf3dds)$stage, levels=c("nashf0f1", "nashf3"))
colData(f0f1vsf4dds)$stage <- factor(colData(f0f1vsf4dds)$stage, levels=c("nashf0f1", "nashf4"))




f0f1vsf2Deseq <- DESeq(f0f1vsf2dds)

f0f1vsf3Deseq <- DESeq(f0f1vsf3dds)

f0f1vsf4Deseq <- DESeq(f0f1vsf4dds)



f0f1vsf2results <- results(f0f1vsf2Deseq)
f0f1vsf3results <- results(f0f1vsf3Deseq)
f0f1vsf4results <- results(f0f1vsf4Deseq)


f0f1vsf2results <- data.frame(f0f1vsf2results)
f0f1vsf3results <- data.frame(f0f1vsf3results)
f0f1vsf4results <- data.frame(f0f1vsf4results)



f0f1vsf2results$gene <- row.names(f0f1vsf2results) 
f0f1vsf3results$gene <- row.names(f0f1vsf3results)
f0f1vsf4results$gene <- row.names(f0f1vsf4results)


                                        #saved original and moving on to make cutoffs


f0f1vsf2results <- f0f1vsf2results[!is.na(f0f1vsf2results$padj),] #18401 initial, 15939 after
f0f1vsf3results <- f0f1vsf3results[!is.na(f0f1vsf3results$padj),] #18401 initial, 16340 after
f0f1vsf4results <- f0f1vsf4results[!is.na(f0f1vsf4results$padj),] #18401 initial, 16480 after


                                        #p<0.001

f0f1vsf2results <- subset(f0f1vsf2results, padj<0.001) #7 DEG vs 0
f0f1vsf3results <- subset(f0f1vsf3results, padj<0.001) #482 DEG vs 434
f0f1vsf4results <- subset(f0f1vsf4results, padj<0.001) #1175 DEG vs 1194

                                        #p<0.01

f0f1vsf2results <- subset(f0f1vsf2results, padj<0.01) #25 DEG vs 0
f0f1vsf3results <- subset(f0f1vsf3results, padj<0.01) #924 DEG vs 434
f0f1vsf4results <- subset(f0f1vsf4results, padj<0.01) #2167 DEG vs 1194


                                        #p<0.005

f0f1vsf2results <- subset(f0f1vsf2results, padj<0.005) #11 DEG vs 0
f0f1vsf3results <- subset(f0f1vsf3results, padj<0.005) #756 DEG vs 434
f0f1vsf4results <- subset(f0f1vsf4results, padj<0.005) #1769 DEG vs 1194

                                        #extracting genes with p<0.001 cutoff

f0f1vsf2DEG <- as.data.frame(f0f1vsf2results[,"gene"], header=T)
f0f1vsf3DEG <- as.data.frame(f0f1vsf3results[,"gene"], header=T)
f0f1vsf4DEG <- as.data.frame(f0f1vsf4results[,"gene"], header=T)


names(f0f1vsf2DEG)[1] <- "NASH-F0-F1vsNASH-F2_DEGs" 
names(f0f1vsf3DEG)[1]  <- "NASH-F0-F1vsNASH-F3_DEGs"
names(f0f1vsf4DEG)[1] <- "NASH-F0-F1vsNASH-F4_DEGs"


write.table(f0f1vsf2DEG, "f0f1vsf2DEG0,005",  sep="\t", col.names=T, row.names=T, quote=F)

naflf0f1results <- read.delim("naflf0f1resultsoriginal", sep="\t", header=T)
naflf2results <- read.delim("naflf2resultsoriginal", sep="\t", header=T)
naflf3results <- read.delim("naflf3resultsoriginal", sep="\t", header=T)
naflf4results <- read.delim("naflf4resultsoriginal", sep="\t", header=T)

f0f1vsf2results <- read.delim("f0f1vsf2resultsoriginal", sep="\t", header=T)
f0f1vsf3results <- read.delim("f0f1vsf3resultsoriginal", sep="\t", header=T)
f0f1vsf4results <- read.delim("f0f1vsf4resultsoriginal", sep="\t", header=T)

write.table(, "",  sep="\t", col.names=T, row.names=T, quote=F)

write.table(, "",  sep="\t", col.names=T, row.names=T, quote=F)



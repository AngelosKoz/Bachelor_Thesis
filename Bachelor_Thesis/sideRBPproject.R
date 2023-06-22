set.seed(1234)

setwd('/home/aggelos/Desktop/Rotation/Files')


protein1 <- read.delim("patients-proteindf", sep="\t", header=T)
protein2 <- as.data.frame(t(protein1))
prot3 <- read.delim("protcod3", sep="\t", header=T)

desconVnafl <- read.delim("controlnaflresultsoriginal", sep="\t", header=T)
desconVf0f1 <- read.delim("controlf0f1resultsoriginal", sep="\t",header=T)


head(desconVnafl)

conVnafl <- protein1[,-c(63:217)]
conVf0f1 <- protein1[,-c(12:62,97:217)]

rbplist <- read.delim("rbplistID", sep="\t",header=T)
rbplistID <- as.data.frame(rbplist[,2])
names(rbplistID) <- "id"
rbplistGENE <- as.data.frame(rbplist[,1])
names(rbplistGENE) <- "gene"



rbp_prot1 <- merge(conVnafl, rbplistID, by="id") #2180
rbp_prot3 <- merge(conVnafl, rbplistGENE, by="gene") #2170

duplicated(rbp_prot1$gene)

rbp_prot2 <- merge(conVf0f1, rbplistID, by="id") #2180
rbp_prot4 <- merge(conVf0f1, rbplistGENE, by="gene") #2170


##See why different ID/GENE try to find unique##
colnames(rbp_prot2)
test1 <- as.data.frame(rbp_prot2[,63])
test2 <- as.data.frame(rbp_prot2[,65])
test3 <- as.data.frame(cbind(test1,test2))
names(test3) <- c("x1","x2")
test2 <- merge(rbp_prot1, rbp_prot2, by="id")
####################END#########################

dessrsf1 <- desconVnafl[grep("SRSF", desconVnafl$id),]
srsf <- conVnafl[grep("SRS", conVnafl$id),]
srsf10 <- conVnafl[grep("SRSF10",conVnafl$id),]
tra <- conVnafl[grep("TRA2",conVnafl$id),]
tra <- tra[-1,]
hnrnp <- conVnafl[grep("HNRNP",conVnafl$id),]
hnrnp
SR <- rbind(srsf,tra)
SR
igf1 <- desconVnafl[grep("IGF1", desconVnafl$id),]
igf1


                                        #SRSF boxplot / heatmap


library(ggplot2)

srsf10 <- srsf10[,-c(1,63,64)]
srsf10t <- as.data.frame(t(srsf10))


srsf10con <- as.data.frame(srsf10t[c(1:10),])
names(srsf10con) <- "counts"
srsf10con$stage <- "Control"

srsf10nafl <- as.data.frame(srsf10t[c(11:61),])
names(srsf10nafl) <- "counts"
srsf10nafl$stage <- "NAFL"

srsf10box <- rbind(srsf10con, srsf10nafl)
srsf10box

ggplot(srsf10box, aes(x=stage, y=counts))+
    geom_boxplot()


                                        #RBP DESeq2

library(DESeq2)


controlstage <- read.delim("controlstage", sep="\t", header=T)
naflstage <- read.delim("naflstage", sep="\t", header=T)
nashf0f1stage <- read.delim("nashf0f1stage", sep="\t", header=T)
controlnafl <- rbind(controlstage, naflstage) #61
controlf0f1 <- rbind(controlstage, nashf0f1stage) #44


head(rbp_prot1) #protein
head(rbp_prot2)
rownames(rbp_prot1) <- rbp_prot1[,2]
rownames(rbp_prot2) <- rbp_prot2[,2]
control_nafl <- rbp_prot1[,-c(1,2,64)]
control_f0f1 <- rbp_prot2[,-c(1,2,47)]


controlnaflframe <- data.frame(stage=c(rep("control", 10), rep('nafl', 51)), replicate=c(1:10, 1:51), type=rep("single",61))
controlf0f1frame <- data.frame(stage=c(rep("control", 10), rep('nashf0f1', 34)), replicate=c(1:10, 1:34), type=rep("single",44))
row.names(controlnaflframe) <- controlnafl[,1]
row.names(controlf0f1frame) <- controlf0f1[,1]


controlnaflCounts <- control_nafl[,match(rownames(controlnaflframe), colnames(control_nafl))]
controlf0f1Counts <- control_f0f1[,match(rownames(controlf0f1frame), colnames(control_f0f1))]


controlnaflCol <- controlnaflframe[,c("stage","type")]
controlf0f1Col <- controlf0f1frame[,c("stage","type")]

controlnafldds <- DESeqDataSetFromMatrix(countData=controlnaflCounts, colData=controlnaflCol, design = ~ stage)
controlf0f1dds <- DESeqDataSetFromMatrix(countData=controlf0f1Counts, colData=controlf0f1Col, design = ~ stage)

colData(controlnafldds)$stage <- factor(colData(controlnafldds)$stage, levels=c("control", "nafl"))
colData(controlf0f1dds)$stage <- factor(colData(controlf0f1dds)$stage, levels=c("control", "nashf0f1"))

controlnaflDeseq <- DESeq(controlnafldds)
controlf0f1Deseq <- DESeq(controlf0f1dds)




controlnaflresults <- results(controlnaflDeseq)
controlnaflresults <- data.frame(controlnaflresults)
controlf0f1results <- results(controlf0f1Deseq)
controlf0f1results <- data.frame(controlf0f1results)

########
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


x1 <- controlnaflCounts[,62:63]
x2 <- controlf0f1Counts[,45:46]

xx <- cbind(x1,x2)

xx <- xx[,-c(3)]

head(xx)

write.table(xx, "GeometricMeansRBP", sep="\t", quote=F)


controlnaflCounts$gene <- row.names(controlnaflCounts)
controlf0f1Counts$gene <- row.names(controlf0f1Counts)

controlnaflresults$gene <- row.names(controlnaflresults)
controlnaflresults <- merge(controlnaflresults, prot3)
controlnaflresults$log10padj <- -log(controlnaflresults$padj,10)
controlf0f1results$gene <- row.names(controlf0f1results)
controlf0f1results <- merge(controlf0f1results, prot3)
controlf0f1results$log10padj <- -log(controlf0f1results$padj,10)

controlnaflresults  <- merge(controlnaflresults, controlnaflCounts[,62:64], by="gene")
controlf0f1results  <- merge(controlf0f1results, controlf0f1Counts[,45:47], by="gene")


head(controlnaflresults)

row.names(controlnaflresults) <- controlnaflresults[,1]
row.names(controlf0f1results) <- controlf0f1results[,1]
###########

dessrsf2 <- controlnaflresults[grep("SRSF10", controlnaflresults$id),]
dessrsf3 <- controlf0f1results[grep("SRSF10", controlf0f1results$id),]


2^(-0.6227903) #find log2
2^(-0.4497529)



write.table(controlnaflresults, "RBP_controlnaflresultsoriginal",  sep="\t", quote = F)
write.table(controlf0f1results, "RBP_controlf0f1resultsoriginal",  sep="\t", quote = F)

controlnaflresults <- controlnaflresults[!is.na(controlnaflresults$padj),] #2180 before, 2163 after 
controlnaflresults <- subset(controlnaflresults, padj<0.005) #2163 before, 821 after
controlnaflresults <- subset(controlnaflresults, padj<0.001) #676


controlf0f1results <- controlf0f1results[!is.na(controlf0f1results$padj),] #2160
controlf0f1results <- subset(controlf0f1results, padj<0.005) #763
controlf0f1results <- subset(controlf0f1results, padj<0.001) #617

nrow(controlf0f1results)


write.table(controlnaflresults, "RBP_controlnaflresults0,001adj",  sep="\t", quote = F)
write.table(controlf0f1results, "RBP_controlf0f1results0,001adj",  sep="\t", quote = F)



                                        #Volcano Plot

library(ggplot2)
library(dplyr)
library(ggrepel)

ctrlVSnafl <- read.delim("RBP_controlnaflresultsoriginal",sep="\t", header=T)
ctrlVSf0f1 <- read.delim("RBP_controlf0f1resultsoriginal",sep="\t", header=T)

connafl <- ctrlVSnafl[!is.na(ctrlVSnafl$padj),] #2163
conf0f1 <- ctrlVSf0f1[!is.na(ctrlVSf0f1$padj),] #2160

unique(connafl$id)
duplicated(connafl$id)
print(connafl$id)


test1 <- connafl[order(connafl$id),] #change parameter
test2 <- rle(test1$id)
test3 <- conf0f1[order(conf0f1$id),] #change parameter
test4 <- rle(test3$id)

test1$id_dupli <- paste0(rep(test2$values, times = test2$lengths), "_",
                         unlist(lapply(test2$lengths, seq_len)))

test1$id_dupli <- gsub("_1", "", test1$id_dupli)

test3$id_dupli <- paste0(rep(test4$values, times = test4$lengths), "_",
                         unlist(lapply(test4$lengths, seq_len)))

test3$id_dupli <- gsub("_1", "", test3$id_dupli)

unique(test1$id_dupli)
duplicated(test1$id_dupli)
print(test1$id_dupli)


connafl <- test1
head(connafl)
conf0f1 <- test3
head(conf0f1)
row.names(connafl)

head(test1)

rownames(connafl) <- connafl[,13]
rownames(conf0f1) <- conf0f1[,13]

connafl[grep("SF3B", connafl$id_dupli),]
tst <- conVnafl[grep("SF3B", connafl$id_dupli),]



wanted  <- as.data.frame(c("SRSF10","ENO1", "ZFP36", "PCBP2", "NOB1","SF3B1"))
names(wanted) <- "id"

commonidlabeling <- merge(wanted, conf0f1, by="id") ##Change DATAFRAME

p01 <- conf0f1 %>% filter(log10padj>=2) ###Change DATAFRAME
p005 <- conf0f1 %>% filter(log10padj>=2.301) ###Change DATAFRAME
p001 <- conf0f1 %>% filter(log10padj>=3) ###Change DATAFRAME


log(1.5,2) #for FC > 1.5
log(0.001,10) #for padj<0.001



ggplot1 <- ggplot(conf0f1, aes(x=log2FoldChange, y=log10padj))+ ###Change DATAFRAME
    scale_x_continuous(breaks = scales::pretty_breaks(n = 15)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    geom_vline(xintercept = 0.5849625, linetype="dashed", color="blue4", alpha=0.5) +
    geom_vline(xintercept = -0.5849625, linetype="dashed", color="blue4", alpha=0.5) +
    geom_hline(yintercept = 2, linetype="dashed", color="blue4", alpha=0.5) +
    geom_hline(yintercept = 2.301, linetype="dashed", color="blue4", alpha=0.5) +
    geom_hline(yintercept = 3, linetype="dashed", color="blue4", alpha=0.5) +
    annotate(x= 0.5849625,y=+Inf,label="1,5 FC",vjust=1, geom="label") +
    annotate(x=- 0.5849625 ,y=+Inf,label="-1,5 FC",vjust=1, geom="label") +
    annotate(x=+Inf ,y=2, label="P=0,01", vjust=1, hjust="inward", geom="label") +
    annotate(x=+Inf ,y=2.301, label="P=0,005", vjust=1, hjust="inward", geom="label") +
    annotate(x=+Inf ,y=3, label="P=0,001", vjust=1, hjust="inward", geom="label") +
    ggtitle("RBP Control vs NASH F0-F1 DESeq2 Results") + ##################################Change NAME
    labs(y="-log10(P-adjusted)", x = "log2(Fold Change)") +
    geom_point(alpha=0.3)+
    geom_point(data=p01, aes(x=log2FoldChange, y=log10padj), color='yellow', size=1)+
    geom_point(data=p005, aes(x=log2FoldChange, y=log10padj), color='green', size=1)+
    geom_point(data=p001, aes(x=log2FoldChange, y=log10padj), color='red', size=1)


ggplot1

#Fix coords

ggplot1 <- ggplot1 + coord_fixed(xlim=c(-5,5), ratio=0.15) ####Fix Coord ylim=c(0,10)

ggplot1

#P<0.01 & abs(FC)>1.5

p <- ggplot1 + geom_text_repel(data=subset(conf0f1, log10padj > 3 & abs(log2FoldChange) > 2), aes(x=log2FoldChange, y=log10padj, label=id), max.overlaps = Inf) ###Change DATAFRAME /  hjust="inward", , vjust="inward",

p

##SPECIFIC
    
p  <- ggplot1 + geom_text_repel(data=subset(connafl, log10padj > 3 & log2FoldChange > 3 | log10padj > 3 & log2FoldChange < -3 |log10padj>3 & abs(log2FoldChange) > 3 | log10padj > 15 & abs(log2FoldChange) > 1.5 ), aes(x=log2FoldChange, y=log10padj, label=id), hjust="inward", , vjust="inward", max.overlaps = Inf) ###Change DATAFRAME

p

##Add the common genes between baseline and stages

p1 <- ggplot1 + geom_text_repel(data=commonidlabeling, aes(x=log2FoldChange, y=log10padj, label=id),  max.overlaps = Inf, colour = "navyblue") ###Change DATAFRAME For LNC / hjust="inward", vjust="outward"

p1








#log(1.5, 2) or log(1/1.5, 2) -> abs(log2(fc)) > 0.5849625



ggsave("RBP-ControlvsNAFL-Scatterplot.pdf",p) ###Change NAME

ggsave("RBP-ControlvsNASHF0-F1-Scatterplot.pdf",p1)



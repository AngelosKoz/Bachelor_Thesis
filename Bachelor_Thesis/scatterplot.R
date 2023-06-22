setwd('/home/aggelos/Desktop/Rotation/Files')

naflf0f1 <- read.delim("naflf0f1resultsoriginal", sep="\t", header=T)
naflf2 <- read.delim("naflf2resultsoriginal", sep="\t", header=T)
naflf3 <-read.delim("naflf3resultsoriginal", sep="\t", header=T)
naflf4 <-read.delim("naflf4resultsoriginal", sep="\t", header=T)

lncnaflf0f1 <- read.delim("lnc_naflf0f1resultsoriginal", sep="\t", header=T) #0 DEG
lncnaflf2 <- read.delim("lnc_naflf2resultsoriginal", sep="\t", header=T)
lncnaflf3 <-read.delim("lnc_naflf3resultsoriginal", sep="\t", header=T)
lncnaflf4 <-read.delim("lnc_naflf4resultsoriginal", sep="\t", header=T)

f0f1vsf2 <- read.delim("f0f1vsf2resultsoriginal", sep="\t", header=T)
f0f1vsf3 <- read.delim("f0f1vsf3resultsoriginal", sep="\t", header=T)
f0f1vsf4 <- read.delim("f0f1vsf4resultsoriginal", sep="\t", header=T)
        
lncf0f1vsf2 <- read.delim("lnc_f0f1vsf2resultsoriginal", sep="\t", header=T) #0
lncf0f1vsf3 <- read.delim("lnc_f0f1vsf3resultsoriginal", sep="\t", header=T)
lncf0f1vsf4 <- read.delim("lnc_f0f1vsf4resultsoriginal", sep="\t", header=T)

prot <- read.delim("protcod3", sep="\t", header=F)
lnc <- read.delim("lncrna3", sep="\t", header=F)
names(lnc) <- c("gene", "id", "type")
names(prot) <- c("gene", "id", "type")




                                    #Removing NA


naflf0f1 <- naflf0f1[!is.na(naflf0f1$padj),] 
naflf2 <- naflf2[!is.na(naflf2$padj),] 
naflf3 <- naflf3[!is.na(naflf3$padj),] 
naflf4 <- naflf4[!is.na(naflf4$padj),]

lncnaflf0f1 <- lncnaflf0f1[!is.na(lncnaflf0f1$padj),] 
lncnaflf2 <- lncnaflf2[!is.na(lncnaflf2$padj),] 
lncnaflf3 <- lncnaflf3[!is.na(lncnaflf3$padj),] 
lncnaflf4 <- lncnaflf4[!is.na(lncnaflf4$padj),]

f0f1vsf2 <- f0f1vsf2[!is.na(f0f1vsf2$padj),] 
f0f1vsf3 <- f0f1vsf3[!is.na(f0f1vsf3$padj),] 
f0f1vsf4 <- f0f1vsf4[!is.na(f0f1vsf4$padj),]

lncf0f1vsf2 <- lncf0f1vsf2[!is.na(lncf0f1vsf2$padj),] 
lncf0f1vsf3 <- lncf0f1vsf3[!is.na(lncf0f1vsf3$padj),] 
lncf0f1vsf4 <- lncf0f1vsf4[!is.na(lncf0f1vsf4$padj),]

                                        #Merge original with ID


naflf0f1 <- merge(naflf0f1, prot, by="gene")
naflf2 <-  merge(naflf2, prot, by="gene")
naflf3 <-  merge(naflf3, prot, by="gene")
naflf4 <-  merge(naflf4, prot, by="gene")

f0f1vsf2 <-  merge(f0f1vsf2, prot, by="gene")
f0f1vsf3 <-  merge(f0f1vsf3, prot, by="gene")
f0f1vsf4 <-  merge(f0f1vsf4, prot, by="gene")

lncnaflf0f1 <- merge(lncnaflf0f1, lnc, by="gene")
lncnaflf2 <-  merge(lncnaflf2, lnc, by="gene")
lncnaflf3 <-  merge(lncnaflf3, lnc, by="gene")
lncnaflf4 <-  merge(lncnaflf4, lnc, by="gene")

lncf0f1vsf2 <-  merge(lncf0f1vsf2, lnc, by="gene")
lncf0f1vsf3 <-  merge(lncf0f1vsf3, lnc, by="gene")
lncf0f1vsf4 <-  merge(lncf0f1vsf4, lnc, by="gene")



                                        #Add ID to rownames


rownames(naflf0f1) <- naflf0f1[,8]
rownames(naflf2) <- naflf2[,8]
rownames(naflf3) <- naflf3[,8]
rownames(naflf4) <- naflf4[,8]

rownames(lncnaflf2) <- lncnaflf2[,8]
rownames(lncnaflf3) <- lncnaflf3[,8]
rownames(lncnaflf4) <- lncnaflf4[,8]

rownames(f0f1vsf2) <- f0f1vsf2[,8]
rownames(f0f1vsf3) <- f0f1vsf3[,8]
rownames(f0f1vsf4) <- f0f1vsf4[,8]

rownames(lncf0f1vsf3) <- lncf0f1vsf3[,8]
rownames(lncf0f1vsf4) <- lncf0f1vsf4[,8]


###Adding Transformed Column of padj with -log10(padj)

naflf0f1$log10padj <- -log(naflf0f1$padj,10)
naflf2$log10padj <- -log(naflf2$padj,10)
naflf3$log10padj <- -log(naflf3$padj,10)
naflf4$log10padj <- -log(naflf4$padj,10)

f0f1vsf2$log10padj <- -log(f0f1vsf2$padj,10)
f0f1vsf3$log10padj <- -log(f0f1vsf3$padj,10)
f0f1vsf4$log10padj <- -log(f0f1vsf4$padj,10)

lncnaflf0f1$log10padj <- -log(lncnaflf0f1$padj,10)
lncnaflf2$log10padj <- -log(lncnaflf2$padj,10)
lncnaflf3$log10padj <- -log(lncnaflf3$padj,10)
lncnaflf4$log10padj <- -log(lncnaflf4$padj,10)

lncf0f1vsf2$log10padj <- -log(lncf0f1vsf2$padj,10)
lncf0f1vsf3$log10padj <- -log(lncf0f1vsf3$padj,10)
lncf0f1vsf4$log10padj <- -log(lncf0f1vsf4$padj,10)




#####Scatterplot

library(ggplot2)
library(dplyr)
library(ggrepel)


commonid <- read.delim("id_nafl0,01", sep="\t", header=T)

commonid <- read.delim("id_f0f1vsf2340,01", sep="\t", header=T)

commonid3 <- read.delim("id_intersectnof20,01", sep="\t", header=T)
commonid4 <- read.delim("id_intersectnof20,001", sep="\t", header=T)


lnccommonid2 <- read.delim("lnc_id_f0f1vsf340,01", sep="\t", header=T)

lnccommonid3 <- read.delim("lnc_id_intersect0,01", sep="\t", header=T)

p01 <- naflf3 %>% filter(log10padj>=2) ###Change DATAFRAME
p005 <- naflf3 %>% filter(log10padj>=2.301) ###Change DATAFRAME
p001 <- naflf3 %>% filter(log10padj>=3) ###Change DATAFRAME

commonidlabeling <- merge(commonid, f0f1vsf4) ###Change DATAFRAME

commonidlabeling2 <- merge(commonid2, naflf3) #(0,001)

lnccommonidlabeling <- merge(lnccommonid2, lncf0f1vsf4) #(0,01)

lnccommonidlabeling2 <- merge(lnccommonid3, lncf0f1vsf2) #(0,001)


ggplot1 <- ggplot(naflf3, aes(x=log2FoldChange, y=log10padj))+ ###Change DATAFRAME
    scale_x_continuous(breaks = scales::pretty_breaks(n = 15)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    geom_vline(xintercept = 1.5, linetype="dashed", color="blue4", alpha=0.5) +
    geom_vline(xintercept = -1.5, linetype="dashed", color="blue4", alpha=0.5) +
    geom_hline(yintercept = 2, linetype="dashed", color="blue4", alpha=0.5) +
    geom_hline(yintercept = 2.301, linetype="dashed", color="blue4", alpha=0.5) +
    geom_hline(yintercept = 3, linetype="dashed", color="blue4", alpha=0.5) +
    annotate(x=1.5 ,y=+Inf,label="1,5 FC",vjust=1, geom="label") +
    annotate(x=-1.5 ,y=+Inf,label="-1,5 FC",vjust=1, geom="label") +
    annotate(x=+Inf ,y=2, label="P=0,01", vjust=1, hjust="inward", geom="label") +
    annotate(x=+Inf ,y=2.301, label="P=0,005", vjust=1, hjust="inward", geom="label") +
    annotate(x=+Inf ,y=3, label="P=0,001", vjust=1, hjust="inward", geom="label") +
    ggtitle("NAFLvsNASH-F3 DESeq2 Results") + ##################################Change NAME
    labs(y="-log10(P-adjusted)", x = "log2(Fold Change)") +
    geom_point(alpha=0.3)+
    geom_point(data=p01, aes(x=log2FoldChange, y=log10padj), color='yellow', size=1)+
    geom_point(data=p005, aes(x=log2FoldChange, y=log10padj), color='green', size=1)+
    geom_point(data=p001, aes(x=log2FoldChange, y=log10padj), color='red', size=1)

ggplot1
    

#Fix coords

ggplot2 <- ggplot1 + coord_fixed(xlim=c(-4,5.5), ratio=0.3) ####Fix Coord ylim=c(0,10)

ggplot2
    

#P<0.01 & abs(FC)>1.5

p <- ggplot1 + geom_text_repel(data=subset(f0f1vsf4, log10padj > 3 & abs(log2FoldChange) > 2), aes(x=log2FoldChange, y=log10padj, label=id), max.overlaps = Inf) ###Change DATAFRAME /  hjust="inward", , vjust="inward",

p

##P<0.01 & abs(FC)>1.5 + P < 0.001
    
p  <- ggplot2 + geom_text_repel(data=subset(naflf3, log10padj > 7 & log2FoldChange > 1.5 | log10padj > 3 & log2FoldChange < -1.5), aes(x=log2FoldChange, y=log10padj, label=id), hjust="inward", , vjust="outward", max.overlaps = Inf) ###Change DATAFRAME

p

##Add the common genes between NAFL and Fibrosis

p1 <- ggplot1 + geom_text_repel(data=commonidlabeling, aes(x=log2FoldChange, y=log10padj, label=id),  max.overlaps = Inf, colour = "navyblue") ###Change DATAFRAME For LNC / hjust="inward", vjust="outward",

p1


ggsave("PCG-NAFLvsF3-Scatterplot2.pdf",p) ###Change NAME

ggsave("PCG-NASHF0F1vsF4-Scatterplot 0,01.pdf",p1) ###Change NAME

ggsave("PCG-NASHF0F1vsF4-Scatterplot 0,001.pdf",p1) ###Change NAME 

                                        #Alt Save

dev.copy(pdf, "PCG-NaflvsF0F1-Scatterplot.pdf") ###Change NAME

dev.off()


head(commonidlabeling)



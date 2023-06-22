setwd('/home/aggelos/Desktop/Rotation/Files')

connafl <- read.delim("controlnaflresultsoriginal", sep="\t", header=T)
conf0f1 <- read.delim("controlf0f1resultsoriginal", sep="\t", header=T)
conf2 <- read.delim("controlf2resultsoriginal", sep="\t", header=T)
conf3 <- read.delim("controlf3resultsoriginal", sep="\t", header=T)
conf4 <- read.delim("controlf4resultsoriginal", sep="\t", header=T)


lncconnafl <- read.delim("lnc_controlnaflresultsoriginal", sep="\t", header=T)
lncconf0f1 <- read.delim("lnc_controlf0f1resultsoriginal", sep="\t", header=T)
lncconf2 <- read.delim("lnc_controlf2resultsoriginal", sep="\t", header=T)
lncconf3 <- read.delim("lnc_controlf3resultsoriginal", sep="\t", header=T)
lncconf4 <- read.delim("lnc_controlf4resultsoriginal", sep="\t", header=T)

connafl <- connafl[!is.na(connafl$padj),]
conf0f1 <- conf0f1[!is.na(conf0f1$padj),]
conf2 <- conf2[!is.na(conf2$padj),]
conf3 <- conf3[!is.na(conf3$padj),]
conf4 <- conf4[!is.na(conf4$padj),]

lncconnafl <- lncconnafl[!is.na(lncconnafl$padj),]
lncconf0f1 <-  lncconf0f1[!is.na(lncconf0f1$padj),]
lncconf2 <- lncconf2[!is.na(lncconf2$padj),]
lncconf3 <- lncconf3[!is.na(lncconf3$padj),]
lncconf4 <- lncconf4[!is.na(lncconf4$padj),]


rownames(connafl) <- connafl[,8]
rownames(conf0f1) <- conf0f1[,8]
rownames(conf2) <- conf2[,8]
rownames(conf3) <- conf3[,8]
rownames(conf4) <- conf4[,8]

rownames(lncconnafl) <- lncconnafl[,8] ##Genes GLIDR/ LINC00864/ TRAF3IP2-AS1 non unique
rownames(lncconf0f1) <- lncconf0f1[,8] ##Genes GLIDR/ LINC00864/ TRAF3IP2-AS1 non unique
rownames(lncconf2) <- lncconf2[,8] ##Genes GLIDR/ LINC00864/ TRAF3IP2-AS1 non unique
rownames(lncconf3) <- lncconf3[,8] ##Genes GLIDR/ LINC00864/ TRAF3IP2-AS1 non unique
rownames(lncconf4) <- lncconf4[,8] ##Genes LINC00864/ TRAF3IP2-AS1 non unique

#####Scatterplot

library(ggplot2)
library(dplyr)
library(ggrepel)

commonid <- read.delim("id_conVSf0f3_0,01", sep="\t", header=T)

commonid2 <- read.delim("id_conVSnafl0,01", sep="\t", header=T)

xy <- c("GAS6-AS1", "MIR4435-1HG", "LINC00665")
wanted <- as.data.frame(xy)
names(wanted)[1] <- "id"

commonidlabeling <- merge(commonid, conf3) ###Change DATAFRAME


commonidlabeling2 <- merge(commonid2, conf4) ###Change DATAFRAME


p01 <- conf3 %>% filter(log10padj>=2) ###Change DATAFRAME
p005 <- conf3 %>% filter(log10padj>=2.301) ###Change DATAFRAME
p001 <- conf3 %>% filter(log10padj>=3) ###Change DATAFRAME

ggplot1 <- ggplot(conf3, aes(x=log2FoldChange, y=log10padj))+ ###Change DATAFRAME
    scale_x_continuous(breaks = scales::pretty_breaks(n = 15)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    geom_vline(xintercept = 1.5, linetype="dashed", color="blue4", alpha=0.5) + #xinter = 1.5
    geom_vline(xintercept = -1.5, linetype="dashed", color="blue4", alpha=0.5) + #xinter = 1.5
    geom_hline(yintercept = 2, linetype="dashed", color="blue4", alpha=0.5) +
    geom_hline(yintercept = 2.301, linetype="dashed", color="blue4", alpha=0.5) +
    geom_hline(yintercept = 3, linetype="dashed", color="blue4", alpha=0.5) +
    annotate(x=1.5 ,y=+Inf,label="1.5 FC",vjust=1, geom="label") + #FC 1.5
    annotate(x=-1.5 ,y=+Inf,label="-1.5 FC",vjust=1, geom="label") + #FC 1.5
    annotate(x=+Inf ,y=2, label="P=0,01", vjust=1, hjust="inward", geom="label") +
    annotate(x=+Inf ,y=2.301, label="P=0,005", vjust=1, hjust="inward", geom="label") +
    annotate(x=+Inf ,y=3, label="P=0,001", vjust=1, hjust="inward", geom="label") +
    ggtitle("ControlvsF3 DESeq2 Results") + ##################################Change NAME
    labs(y="-log10(P-adjusted)", x = "log2(Fold Change)") +
    geom_point(alpha=0.01)+ #geom alpha 0.3 / 0.01
    geom_point(data=p01, aes(x=log2FoldChange, y=log10padj), color='yellow', size=1)+
    geom_point(data=p005, aes(x=log2FoldChange, y=log10padj), color='green', size=1)+
    geom_point(data=p001, aes(x=log2FoldChange, y=log10padj), color='red', size=1)

ggplot1

#Fix coords

ggplot1 <- ggplot1 + coord_fixed(xlim=c(-8.5,8.5), ratio=0.15) ####Fix Coord ylim=c(0,10)

ggplot1

#P<0.01 & abs(FC)>1.5

p <- ggplot1 + geom_text_repel(data=subset(lncconf0f1, log10padj > 3 & abs(log2FoldChange) > 2), aes(x=log2FoldChange, y=log10padj, label=id), max.overlaps = Inf) ###Change DATAFRAME /  hjust="inward", , vjust="inward",

p

##SPECIFIC
    
p  <- ggplot1 + geom_text_repel(data=subset(lncconnafl, log10padj > 3 & log2FoldChange > 3 | log10padj > 3 & log2FoldChange < -3 |log10padj>3 & abs(log2FoldChange) > 3 | log10padj > 15 & abs(log2FoldChange) > 1.5 ), aes(x=log2FoldChange, y=log10padj, label=id), hjust="inward", , vjust="outward", max.overlaps = Inf) ###Change DATAFRAME

p

##Add the common genes between baseline and stages

p1 <- ggplot1 + geom_text_repel(data=commonidlabeling, aes(x=log2FoldChange, y=log10padj, label=id),  max.overlaps = Inf, colour = "navyblue") ###Change DATAFRAME For LNC / hjust="inward", vjust="outward",

p1




ggsave("LNC-ControlvsF4-Scatterplot.pdf",p) ###Change NAME

ggsave("PCG-ControlvsF3-Scatterplot-f0f30,01.pdf",p1)

ggsave("PCG-ControlvsF4-Scatterplot-nafl0,01.pdf",p1)

setwd('/home/aggelos/Desktop/Rotation/Files')

set.seed(1234)

library(pls)
library(pheatmap)
library(dplyr)

lnc <- read.delim("lncrna", sep="\t", header=T)

y1 <- read.delim("lnc_controlnaflresults0,01adj", sep="\t", header = T)
y2 <- read.delim("lnc_controlf0f1results0,01adj", sep="\t", header = T)
y3 <- read.delim("lnc_controlf2results0,01adj", sep="\t", header = T)
y4 <- read.delim("lnc_controlf3results0,01adj", sep="\t", header = T)
y5 <- read.delim("lnc_controlf4results0,01adj", sep="\t", header = T)

head(y1)

nrow(y5)

y1 <- subset(y1, abs(log2FoldChange) > 0.5849625)
y2 <- subset(y2, abs(log2FoldChange) > 1.5)
y3 <- subset(y3, abs(log2FoldChange) > 1.5)
y4 <- subset(y4, abs(log2FoldChange) > 1.5)
y5 <- subset(y5, abs(log2FoldChange) > 1.5)

y1 <- y1[, c(1,8)]
y2 <- y2[, c(1,8)]
y3 <- y3[, c(1,8)]
y4 <- y4[, c(1,8)]
y5 <- y5[, c(1,8)]

degs0 <- rbind(y1,y2,y3,y4,y5)
degs <- distinct(degs0) #alternate <- distinct(degs0, gene, .keep_all=T)



##Patients Heatmap

degs <- merge(degs, lnc)
row.names(degs) <- degs[,2]


controlstage <- read.delim("controlstage", sep="\t", header=T)
naflstage <- read.delim("naflstage", sep="\t", header=T)
nashf0f1stage <- read.delim("nashf0f1stage", sep="\t", header=T)
nashf2stage <- read.delim("nashf2stage", sep="\t", header=T)
nashf3stage <- read.delim("nashf3stage", sep="\t", header=T)
nashf4stage <- read.delim("nashf4stage", sep="\t", header=T)

control <- degs[, as.character(names(degs)) %in% as.character(controlstage$patient)]
nafl <- degs[, as.character(names(degs)) %in% as.character(naflstage$patient)]
f0f1 <- degs[, as.character(names(degs)) %in% as.character(nashf0f1stage$patient)]
f2 <- degs[, as.character(names(degs)) %in% as.character(nashf2stage$patient)]
f3 <- degs[, as.character(names(degs)) %in% as.character(nashf3stage$patient)]
f4 <- degs[, as.character(names(degs)) %in% as.character(nashf4stage$patient)]

df <- cbind(control, nafl, f0f1, f2, f3, f4)
df1 <- as.matrix(df)
df1 <- stdize(df1,center = T, scale = T) #center/scale = T?
df1 <- as.data.frame(df1)

pheatmap(df1, main = "Control Baseline Heatmap", cluster_rows=T, cluster_col=F)



dftest <- df
dftest[dftest==0] <- 1 ##make all 0 into 1
dftest1 <- log(dftest, 10)

dftest2 <- as.matrix(dftest1)
dftest2 <- stdize(dftest2, scale=T, center = T)

dftest2 <- as.data.frame(dftest2)
max(dftest2)
min(dftest2)

pheatmap(dftest2, main = "Control Baseline Heatmap", cluster_rows=T, cluster_col=F)

pheatmap(dftest2, main = "Control Baseline Heatmap", cluster_rows=T, cluster_col=F,plot=F)

dev.copy(pdf, "test1.pdf")
dev.off()


##try ggplot heatmap
##pdf- illustrator
##consider log scaling and/or standardize


##Geometric Mean Heatmap


geomean <- read.delim("GeometricMeansLNC", sep="\t", header=T)
geomean$gene <- row.names(geomean)

geodegs <- merge(degs, geomean, by = "gene")
row.names(geodegs) <- geodegs[,2]


geodegs1 <- geodegs[,-c(1,2)]

geodegs2 <- log(geodegs1, 10)

geodegs2 <- as.matrix(geodegs2)
geodegs2 <- stdize(geodegs2, scale=T, center = T)
geodegs2 <- as.data.frame(geodegs2)

max(geodegs2)
min(geodegs2)

nrow(geodegs2)

pheatmap(geodegs2, main = "Control Baseline Heatmap", cluster_rows=T, cluster_col=F)

dev.copy(pdf, "PheatGeo(abs(log2f)c>1.5).pdf")
dev.off()

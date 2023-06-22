setwd('/home/aggelos/Desktop/Rotation/Files')

library(pls)
library(ROCR)
library(ggplot2)
library(ggbiplot)

library(pca3d)

library(rgl)

library(dplyr)

controlstage <- read.delim("controlstage", sep="\t", header=T)
naflstage <- read.delim("naflstage", sep="\t", header=T)
f0f1stage <- read.delim("nashf0f1stage", sep="\t", header=T)
f2stage <- read.delim("nashf2stage", sep="\t", header=T)
f3stage <- read.delim("nashf3stage", sep="\t", header=T)
f4stage <- read.delim("nashf4stage", sep="\t", header=T)

degsf2 <- read.delim("naflf2results0,01adj", sep="\t", header=T)
degsf3 <- read.delim("naflf3results0,01adj", sep="\t", header=T)
degsf4 <- read.delim("naflf4results0,01adj", sep="\t", header=T)

conpat <- read.delim("conpat", sep="\t", header=T)
rownames(conpat) <- conpat[,1]
control <- conpat[,1:11]

#find the corresponding patients in the row data

nafl <- conpat[, as.character(names(conpat)) %in% as.character(naflstage$patient)]
nafl <- tibble::rownames_to_column(nafl, "gene")
rownames(nafl) <- nafl[,1]

nashf0f1 <- conpat[, as.character(names(conpat)) %in% as.character(f0f1stage$patient)]
nashf0f1 <- tibble::rownames_to_column(nashf0f1, "gene")
rownames(nashf0f1) <- nashf0f1[,1]

nashf2 <- conpat[, as.character(names(conpat)) %in% as.character(f2stage$patient)]
nashf2 <- tibble::rownames_to_column(nashf2, "gene")
rownames(nashf2) <- nashf2[,1]

nashf3 <- conpat[, as.character(names(conpat)) %in% as.character(f3stage$patient)]
nashf3 <- tibble::rownames_to_column(nashf3, "gene")
rownames(nashf3) <- nashf3[,1]

nashf4 <- conpat[, as.character(names(conpat)) %in% as.character(f4stage$patient)]
nashf4 <- tibble::rownames_to_column(nashf4, "gene")
rownames(nashf4) <- nashf4[,1]


##control_transpose = t(control_original)

controlx <- control[,-1]
controlx <- as.data.frame(t(controlx)) ##transform columns in to rows

naflx <- nafl[,-1]
naflx <- as.data.frame(t(naflx))

nashf0f1x <- nashf0f1[,-1]
nashf0f1x <- as.data.frame(t(nashf0f1x))

nashf2x <- nashf2[,-1]
nashf2x <- as.data.frame(t(nashf2x))

nashf3x <- nashf3[,-1]
nashf3x <- as.data.frame(t(nashf3x))

nashf4x <- nashf4[,-1]
nashf4x <- as.data.frame(t(nashf4x))


##Filter the specific DEG set



controlmod <- controlx[, as.character(names(controlx)) %in% as.character(degsf3$gene)] ##Change degs for comparisons
naflmod <- naflx[, as.character(names(naflx)) %in% as.character(degsf3$gene)]
nashf0f1mod <- nashf0f1x[, as.character(names(nashf0f1x)) %in% as.character(degsf3$gene)]
nashf2mod <- nashf2x[, as.character(names(nashf2x)) %in% as.character(degsf3$gene)]
nashf3mod <- nashf3x[, as.character(names(nashf3x)) %in% as.character(degsf3$gene)]
nashf4mod <- nashf4x[, as.character(names(nashf4x)) %in% as.character(degsf3$gene)]


#take random columns ##START PCA 

randomnafl <- naflmod[sample(nrow(naflmod), 10), ]
randomf0f1 <- nashf0f1mod[sample(nrow(nashf0f1mod), 10), ]
randomf2 <- nashf2mod[sample(nrow(nashf2mod), 10), ]
randomf3 <- nashf3mod[sample(nrow(nashf3mod), 10), ]
randomf4 <- nashf4mod[sample(nrow(nashf4mod), 10), ]

controlmod$stage  <- "Control"
randomnafl$stage <- "NAFL"
randomf0f1$stage  <- "NASHF0F1"
randomf2$stage  <- "NASHF2"
randomf3$stage  <- "NASHF3"
randomf4$stage  <- "NASHF4"

combined1 <- rbind(controlmod, randomf0f1, randomf2)

pca <- prcomp(combined1[, -1792],  scale = TRUE)
plot <- ggbiplot(pca, groups=factor (combined1$stage), alpha=1, varname.size=0.2, ellipse = TRUE, circle = TRUE, labels =row.names(combined1)) +theme_bw()
plot


plot1 <- ggbiplot(pca, groups=factor (combined1$stage), alpha=1, varname.size=0.2, ellipse = TRUE, circle = TRUE, obs.scale = 1, var.scale = 1)+theme_bw()
plot1

ggsave("PCA controlvsf2.pdf",plot)


##Reminder to filter the up-regulated genes and try PCA again. filter from the deseq results
##Take random 20-30 patients no matter the stage

vstTctrlf2 <- read.delim("transvstTRUE", sep="\t", header=T)
vstFctrlf2 <- read.delim("transvstFALSE", sep="\t", header=T)

vstTctrlf2 <- as.data.frame(t(vstTctrlf2))
vstFctrlf2 <- as.data.frame(t(vstFctrlf2))
nrow(vstFctrlf2)
ncol(vstFctrlf2)

###vstTnaflf2mod <- vstTnaflf2[, as.character(names(vstTnaflf2)) %in% as.character(degsf$gene)]
###vstFnaflf2mod <- vstFnaflf2[, as.character(names(vstFnaflf2)) %in% as.character(degsf$gene)]

vstFctrlf2 <- vstFctrlf2[,apply(vstFctrlf2, 2, var, na.rm=TRUE) != 0] ##remove 0 variance columns -> 18401 before/ 17940 after
ncol(vstFctrlf2)

tst <- vstFctrlf2[(1:10),]
tst1 <- vstFctrlf2[c(11:62),]
tst$stage <- "Control"
tst1$stage <- "NASH-F2"
tst <- as.data.frame(tst[,"stage"], header=T)
tst1 <- as.data.frame(tst1[,"stage"], header=T)
names(tst)[1] <- "stage"
names(tst1)[1] <- "stage"

tst2 <- rbind(tst,tst1)
vstFctrlf2 <- cbind(tst2, vstFctrlf2)
ncol(vstFctrlf2)


write.table(vstFctrlf2, "PCActrlF2data", quote=F, sep="\t")

vstFctrlf2 <- read.delim("PCActrlF2data", sep="\t", header=T)

head(vstFctrlf2, 1)

pca <- prcomp(vstFctrlf2[, -1],  scale = FALSE)

summary(pca)

plot <- ggbiplot(pca, groups=factor(vstFctrlf2$stage), alpha=1, varname.size=0.2, ellipse = TRUE, circle = TRUE, labels =row.names(vstFctrlf2)) +theme_bw()
plot

plot2 <- ggbiplot(pca,choices = c(1:2), groups=factor (vstFctrlf2$stage), alpha=1, varname.size=0.2, ellipse = TRUE, circle = TRUE, labels =row.names(vstFctrlf2)) +theme_bw()
plot2

plot3d(pca$scores[,1:3], col="red", type="p")



##Antwnis

library("scatterplot3d")

source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')

### pca is your PCA values , colors_covid is a vector with an equal number of color names corresponding to your sample. Example : If my samples are Healthy, Healthy,Patient,Super Patient,Patient my color vector can be colors_covid=c("grey","grey","steelblue 4","firebrick 2","steelblue 4")


scatterplot3d(cex.symbols=2.5,x=pca$x[,1],pca$x[,2],pca$x[,3],color=?? ,pch=19,grid=F, font.axis=2,cex.lab=1.2, font.main=2, cex.main=1.5, box=F
,xlab="PC1 ~ 20%", ylab="PC2 ~ 13%", zlab="PC3",main="PCA analysis of DESeq2 normalized counts")
addgrids3d(x=pca$x[,1],pca$x[,2],pca$x[,3], grid = c("xy", "xz", "yz"))
legend("bottom", legend = c("Control","NASH-F2"),
       col =  c("steelblue 4","orange"), pch = 16,inset = -0.18, xpd = TRUE, horiz = TRUE, pt.cex=2,text.font=2)



##Not working
stage <- factor(vstFctrlf2[,1])
summary(stage)

pca3d(pca, group=stage, show.ellipses=TRUE,
ellipse.ci=0.75, show.plane=FALSE)

snapshotPCA3d(file="CtrlvsF2_3D.png")


##Full Antwnis

scatterplot3d(cex.symbols=2.5,x=pca$x[,1],pca$x[,2],pca$x[,3],color=colors_covid,pch=19,grid=F,font.axis=2,cex.lab=1.2,font.main=2,cex.main=1.5,box=F
,xlab="PC1 ~ 37%", ylab="PC2 ~ 12,5%", zlab="PC3 ~ 8,7%",main="PCA analysis of DESeq2 normalized counts")
addgrids3d(x=pca$x[,1],pca$x[,2],pca$x[,3], grid = c("xy", "xz", "yz"))
legend("bottom", legend = c("Control","Covid","FM","Covid-FM"),
       col =  c("steelblue 4", "grey", "yellow","orange"), pch = 16,inset = -0.18, xpd = TRUE, horiz = TRUE, pt.cex=2,text.font=2)

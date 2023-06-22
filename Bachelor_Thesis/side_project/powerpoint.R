setwd('/home/aggelos/Desktop/Rotation/Files')

x1 <- read.delim("lnc_id_nafl0,01", sep="\t", header = T)
x2 <- read.delim("lnc_id_nafl0,005", sep="\t", header = T)
x3 <- read.delim("lnc_id_nafl0,001", sep="\t", header = T)

x4 <- read.delim("lnc_id_f0f1vsf340,01", sep="\t", header = T)
x5 <- read.delim("lnc_id_f0f1vsf340,005", sep="\t", header = T)
x6 <- read.delim("lnc_id_f0f1vsf340,001", sep="\t", header = T)

x7 <- read.delim("lnc_id_intersect0,01", sep="\t", header = T)
x8 <- read.delim("lnc_id_intersect0,005", sep="\t", header = T)
x9 <- read.delim("lnc_id_intersect0,001", sep="\t", header = T)

y1 <- read.delim("id_nafl0,01", sep="\t", header = T)
y2 <- read.delim("id_nafl0,005", sep="\t", header = T)
y3 <- read.delim("id_nafl0,001", sep="\t", header = T)

y4 <- read.delim("id_f0f1vsf340,01", sep="\t", header = T)
y5 <- read.delim("id_f0f1vsf340,005", sep="\t", header = T)
y6 <- read.delim("id_f0f1vsf340,001", sep="\t", header = T)

y7 <- read.delim("id_intersectnof20,01", sep="\t", header = T)
y8 <- read.delim("id_intersectnof20,005", sep="\t", header = T)
y9 <- read.delim("id_intersectnof20,001", sep="\t", header = T)

z4 <- read.delim("id_f0f1vsf2340,01", sep="\t", header = T)
z5 <- read.delim("id_f0f1vsf2340,005", sep="\t", header = T)
z6 <- read.delim("id_f0f1vsf2340,001", sep="\t", header = T)

z7 <- read.delim("id_intersectwithf20,01", sep="\t", header = T)
z8 <- read.delim("id_intersectwithf20,005", sep="\t", header = T)
z9 <- read.delim("id_intersectwithf20,001", sep="\t", header = T)

paper <- read.delim("paperdeg", sep="\t", header=T)

paper1  <- as.data.frame( paper[, 1])

paper1

papercommon <- merge(paper, y4)

papercommon1 <- as.data.frame(papercommon[,1])

papercommon1

paper1

x1 <- x1[2]
x2 <- x2[2]
x3 <- x3[2]

x4 <- x4[2]
x5 <- x5[2]
x6 <- x6[2]

x7 <- x7[2]
x8 <- x8[2]
x9 <- x9[2]

y1 <- y1[2]
y2 <- y2[2]
y3 <- y3[2]

y4 <- y4[2]
y5 <- y5[2]
y6 <- y6[2]

y7 <- y7[2]
y8 <- y8[2]
y9 <- y9[2]

z4 <- z4[2]
z5 <- z5[2]
z6 <- z6[2]

z7 <- z7[2]
z8 <- z8[2]
z9 <- z9[2]

z6 <- z6[order(z6$id),]
z6 <- as.data.frame(z6)


                                        #pheatmap transformations


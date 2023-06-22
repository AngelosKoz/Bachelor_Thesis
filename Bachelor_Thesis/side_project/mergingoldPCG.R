setwd('/home/aggelos/Desktop/Rotation/Files')

prot <- read.delim("protcod3", sep="\t", header=F)
names(prot) <- c("gene", "id", "type")
prot <- prot[,-3]

nafl001 <- read.delim("naflvsf234common0,001", sep="\t", header=F)
nafl001 <- nafl001[-1,]
names(nafl001) <- c("number", "gene")

mernafl001 <- merge(nafl001, prot, by="gene")
mernafl001 <- mernafl001[,-2]


nafl005 <- read.delim("naflvsf234common0,005", sep="\t", header=F)
nafl005 <- nafl005[-1,]
names(nafl005) <- c("number", "gene")
mernafl005 <- merge(nafl005, prot, by="gene")
mernafl005 <- mernafl005[,-2]


nafl01 <- read.delim("naflvsf234common0,01", sep="\t", header=F)
nafl01 <- nafl01[-1,]
names(nafl01) <- c("number", "gene")
mernafl01 <- merge(nafl01, prot, by="gene")
mernafl01 <- mernafl01[,-2]



f0f1vsf234001 <- read.delim("f0f1vsf234common0,001", sep="\t", header=F)
f0f1vsf234001 <- f0f1vsf234001[-1,]
names(f0f1vsf234001) <- c("number", "gene")
merf0f1vsf234001 <- merge(f0f1vsf234001, prot, by="gene")
merf0f1vsf234001 <- merf0f1vsf234001[,-2]


f0f1vsf234005 <- read.delim("f0f1vsf234common0,005", sep="\t", header=F)
f0f1vsf234005 <- f0f1vsf234005[-1,]
names(f0f1vsf234005) <- c("number", "gene")
merf0f1vsf234005 <- merge(f0f1vsf234005, prot, by="gene")
merf0f1vsf234005 <- merf0f1vsf234005[,-2]


f0f1vsf23401 <- read.delim("f0f1vsf234common0,01", sep="\t", header=F)
f0f1vsf23401 <- f0f1vsf23401[-1,]
names(f0f1vsf23401) <- c("number", "gene")
merf0f1vsf23401 <- merge(f0f1vsf23401, prot, by="gene")
merf0f1vsf23401 <- merf0f1vsf23401[,-2]


f0f1vsf34001 <- read.delim("f0f1vsf34common0,001", sep="\t", header=F)
f0f1vsf34001 <- f0f1vsf34001[-1,]
names(f0f1vsf34001) <- c("number", "gene")
merf0f1vsf34001 <- merge(f0f1vsf34001, prot, by="gene")
merf0f1vsf34001 <- merf0f1vsf34001[,-2]


f0f1vsf34005 <- read.delim("f0f1vsf34common0,005", sep="\t", header=F)
f0f1vsf34005 <- f0f1vsf34005[-1,]
names(f0f1vsf34005) <- c("number", "gene")
merf0f1vsf34005 <- merge(f0f1vsf34005, prot, by="gene")
merf0f1vsf34005 <- merf0f1vsf34005[,-2]


f0f1vsf3401 <- read.delim("f0f1vsf34common0,01", sep="\t", header=F)
f0f1vsf3401 <- f0f1vsf3401[-1,]
names(f0f1vsf3401) <- c("number", "gene")
merf0f1vsf3401 <- merge(f0f1vsf3401, prot, by="gene")
merf0f1vsf3401 <- merf0f1vsf3401[,-2]

##common nafl/f0f1 baselines

idnafl001 <- read.delim("id_nafl0,001", sep="\t", header=T)
idnafl005 <- read.delim("id_nafl0,005", sep="\t", header=T)
idnafl01 <- read.delim("id_nafl0,01", sep="\t", header=T)

idf0v234001 <- read.delim("id_f0f1vsf2340,001", sep="\t", header=T)
idf0v234005 <- read.delim("id_f0f1vsf2340,005", sep="\t", header=T)
idf0v23401 <- read.delim("id_f0f1vsf2340,01", sep="\t", header=T)

idf0v34001 <- read.delim("id_f0f1vsf340,001", sep="\t", header=T)
idf0v34005 <- read.delim("id_f0f1vsf340,005", sep="\t", header=T)
idf0v3401 <- read.delim("id_f0f1vsf340,01", sep="\t", header=T)

mernaflf0f2001 <- merge(idnafl001, idf0v234001)
mernaflf0f2001


mernaflf0f2005 <- merge(idnafl005, idf0v234005)
mernaflf0f2005

mernaflf0f201 <- merge(idnafl01, idf0v23401)
mernaflf0f201

mernaflf0f3001 <- merge(idnafl001, idf0v34001)
mernaflf0f3001

mernaflf0f3005 <- merge(idnafl005, idf0v34005)
mernaflf0f3005


mernaflf0f301 <- merge(idnafl01, idf0v3401)
mernaflf0f301

#Checking paper DEG

id <- c("AKR1B10", "ANKRD29", "CCL20", "CFAP221", "CLIC6", "COL1A1", "COL1A2", "DTNA", "DUSP8", "EPB41L4A", "FERMT1", "GDF15", "GDF15", "HSD17B14", "IL32", "TGBL1", "LTBP2", "PDGFA", "PPAPDC1A", "RGS4", "SCTR", "STMN2", "THY1", "TNFRSF12A", "TYMS")
paperDEG <- data.frame(id)

write.table(paperDEG, "paperDEG", sep="\t", col.names=T, row.names=T, quote=F)

check <- merge(prot, paperDEG)
check

merpaper001 <- merge(paperDEG, mernaflf0f3001)
merpaper001

merpaper005 <- merge(paperDEG, mernaflf0f3005)
merpaper005

merpaper01 <- merge(paperDEG, mernaflf0f301)
merpaper01


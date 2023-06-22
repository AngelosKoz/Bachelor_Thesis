setwd('/home/aggelos/Desktop/Rotation/Files')

x1 <- read.delim("naflf2DEG0,001", header=T, sep="\t")
x2 <- read.delim("naflf2DEG0,005", header=T, sep="\t")
x3 <- read.delim("naflf2DEG0,01", header=T, sep="\t")

x4 <- read.delim("naflf3DEG0,001", header=T, sep="\t")
x5 <- read.delim("naflf3DEG0,005", header=T, sep="\t")
x6 <- read.delim("naflf3DEG0,01", header=T, sep="\t")

x7 <- read.delim("naflf4DEG0,001", header=T, sep="\t")
x8 <- read.delim("naflf4DEG0,005", header=T, sep="\t")
x9 <- read.delim("naflf4DEG0,01", header=T, sep="\t")

x10 <- read.delim("f0f1vsf2DEG0,001", header=T, sep="\t")
x11 <- read.delim("f0f1vsf2DEG0,005", header=T, sep="\t")
x12 <- read.delim("f0f1vsf2DEG0,01", header=T, sep="\t")

x13 <- read.delim("f0f1vsf4DEG0,001", header=T, sep="\t")
x14 <- read.delim("f0f1vsf4DEG0,005", header=T, sep="\t")
x15 <- read.delim("f0f1vsf4DEG0,01", header=T, sep="\t")

x16 <- read.delim("f0f1vsf3DEG0,001", header=T, sep="\t")
x17 <- read.delim("f0f1vsf3DEG0,005", header=T, sep="\t")
x18 <- read.delim("f0f1vsf3DEG0,01", header=T, sep="\t")

x19 <- read.delim("lnc_naflf2DEG0,005", header=T, sep="\t")
x20 <- read.delim("lnc_naflf2DEG0,01", header=T, sep="\t")
x21 <- read.delim("lnc_naflf2DEG0,001", header=T, sep="\t")

x22 <- read.delim("lnc_naflf3DEG0,005", header=T, sep="\t")
x23 <- read.delim("lnc_naflf3DEG0,01", header=T, sep="\t")
x24 <- read.delim("lnc_naflf3DEG0,001", header=T, sep="\t")

x25 <- read.delim("lnc_naflf4DEG0,001", header=T, sep="\t")
x26 <- read.delim("lnc_naflf4DEG0,005", header=T, sep="\t")
x27 <- read.delim("lnc_naflf4DEG0,01", header=T, sep="\t")

x28 <- read.delim("lnc_f0f1vsf3DEG0,001", header=T, sep="\t")
x29 <- read.delim("lnc_f0f1vsf3DEG0,005", header=T, sep="\t")
x30 <- read.delim("lnc_f0f1vsf3DEG0,01", header=T, sep="\t")

x31 <- read.delim("lnc_f0f1vsf4DEG0,001", header=T, sep="\t")
x32 <- read.delim("lnc_f0f1vsf4DEG0,005", header=T, sep="\t")
x33 <- read.delim("lnc_f0f1vsf4DEG0,01", header=T, sep="\t")

lnc <- read.delim("lncrna3", header=T, sep="\t")
prot <- read.delim("protcod3", header=T, sep="\t")

x1 <- merge(prot, x1)
x2 <- merge(prot, x2)
x3 <- merge(prot, x3)

nrow(x3)

x4 <- merge(prot, x4, by="gene")
x5 <- merge(prot, x5, by="gene")
x6 <- merge(prot, x6, by="gene")

x7 <- merge(prot, x7)
x8 <- merge(prot, x8)
x9 <- merge(prot, x9)

x10 <- merge(prot, x10)
x11 <- merge(prot, x11)
x12 <- merge(prot, x12)

x13 <- merge(prot, x13)
x14 <- merge(prot, x14)
x15 <- merge(prot, x15)

x16 <- merge(prot, x16)
x17 <- merge(prot, x17)
x18 <- merge(prot, x18)

x19 <- merge(lnc, x19)
x20 <- merge(lnc, x20)
x21 <- merge(lnc, x21)

x22 <- merge(lnc, x22)
x23 <- merge(lnc, x23)
x24 <- merge(lnc, x24)

x25 <- merge(lnc, x25)
x26 <- merge(lnc, x26)
x27 <- merge(lnc, x27)

x28 <- merge(lnc, x28)
x29 <- merge(lnc, x29)
x30 <- merge(lnc, x30)

x31 <- merge(lnc, x31)
x32 <- merge(lnc, x32)
x33 <- merge(lnc, x33)




    
write.table(x1, "naflf2DEG0,001", sep="/t", quote=F)
write.table(x2, "naflf2DEG0,005", sep="/t", quote=F)
write.table(x3, "naflf2DEG0,01", sep="/t", quote=F)

write.table(x4, "naflf3DEG0,001", sep="/t", quote=F)
write.table(x5, "naflf3DEG0,005", sep="/t", quote=F)
write.table(x6, "naflf3DEG0,01", sep="/t", quote=F)

write.table(x7, "naflf4DEG0,001", sep="/t", quote=F)
write.table(x8, "naflf4DEG0,005", sep="/t", quote=F)
write.table(x9, "naflf4DEG0,01", sep="/t", quote=F)

write.table(x10, "f0f1vsf2DEG0,001", sep="/t", quote=F)
write.table(x11, "f0f1vsf2DEG0,005", sep="/t", quote=F)
write.table(x12, "f0f1vsf2DEG0,01", sep="/t", quote=F)

write.table(x13, "f0f1vsf4DEG0,001", sep="/t", quote=F)
write.table(x14, "f0f1vsf4DEG0,005", sep="/t", quote=F)
write.table(x15, "f0f1vsf4DEG0,01", sep="/t", quote=F)

write.table(x16, "f0f1vsf3DEG0,001", sep="/t", quote=F)
write.table(x17, "f0f1vsf3DEG0,005", sep="/t", quote=F)
write.table(x18, "f0f1vsf3DEG0,01", sep="/t", quote=F)

write.table(x19, "lnc_naflf2DEG0,005", sep="/t", quote=F)
write.table(x20, "lnc_naflf2DEG0,01", sep="/t", quote=F)
write.table(x21, "lnc_naflf2DEG0,001", sep="/t", quote=F)
                 
write.table(x22, "lnc_naflf3DEG0,005", sep="/t", quote=F)
write.table(x23, "lnc_naflf3DEG0,01", sep="/t", quote=F)
write.table(x24, "lnc_naflf3DEG0,001", sep="/t", quote=F)

write.table(x25, "lnc_naflf4DEG0,001", sep="/t", quote=F)
write.table(x26, "lnc_naflf4DEG0,005", sep="/t", quote=F)
write.table(x27, "lnc_naflf4DEG0,01", sep="/t", quote=F)

write.table(x28, "lnc_f0f1vsf3DEG0,001", sep="/t", quote=F)
write.table(x29, "lnc_f0f1vsf3DEG0,005", sep="/t", quote=F)
write.table(x30, "lnc_f0f1vsf3DEG0,01", sep="/t", quote=F)

write.table(x31, "lnc_f0f1vsf4DEG0,001", sep="/t", quote=F)
write.table(x32, "lnc_f0f1vsf4DEG0,005", sep="/t", quote=F)
write.table(x33, "lnc_f0f1vsf4DEG0,01", sep="/t", quote=F)

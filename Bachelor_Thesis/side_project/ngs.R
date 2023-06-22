x <- numeric(55)
x[1]  <- 0
x[2] <- 1
for(i in 3:55) x[i] <- x[i-2]+x[i-1]

answer <- (x[55]-x[45])*x[5]
answer

#####



sequence <-

most <- function(x) {
  tab <- table(strsplit(x, '')[[1]])
  names(tab)[tab == max(tab)]}

least <- function(x) {
  tab <- table(strsplit(x, '')[[1]])
  names(tab)[tab == min(tab)]}

most(sequence)
lengths(regmatches(sequence, gregexpr("S", sequence)))
least(sequence)
lengths(regmatches(sequence, gregexpr("W", sequence)))

#####

reference <- ("ATGACCCCAATACGCAAAACTAACCCCCTAATAAAATTAATTAACCACTCATTCATCGACCTCCCCACCCCATCCAACATCTCCGCATGATGAAACTTCGGCTCACTCCTTGGCGCCTGCCTGATCCTCCAAATCACCACAGGACTATTCCTAGCCATGCACTACTCACCAGACGCCTCAACCGCCTTTTCATCAATCGCCCACATCACTCGAGACGTAAATTATGGCTGAATCATCCGCTACCTTCACGCCAATGGCGCCTCAATATTCTTTATCTGCCTCTTCCTACACATCGGGCGAGGCCTATATTACGGATCATTTCTCTACTCAGAAACCTGAAACATCGGCATTATCCTCCTGCTTGCAACTATAGCAACAGCCTTCATAGGCTATGTCCTCCCGTGAGGCCAAATATCATTCTGAGGGGCCACAGTAATTACAAACTTACTATCCGCCATCCCATACATTGGGACAGACCTAGTTCAATGAATCTGAGGAGGCTACTCAGTAGACAGTCCCACCCTCACACGATTCTTTACCTTTCACTTCATCTTGCCCTTCATTATTGCAGCCCTAGCAACACTCCACCTCCTATTCTTGCACGAAACGGGATCAAACAACCCCCTAGGAATCACCTCCCATTCCGATAAAATCACCTTCCACCCTTACTACACAATCAAAGACGCCCTCGGCTTACTTCTCTTCCTTCTCTCCTTAATGACATTAACACTATTCTCACCAGACCTCCTAGGCGACCCAGACAATTATACCCTAGCCAACCCCTTAAACACCCCTCCCCACATCAAGCCCGAATGATATTTCCTATTCGCCTACACAATTCTCCGATCCGTCCCTAACAAACTAGGAGGCGTCCTTGCCCTATTACTATCCATCCTCATCCTAGCAATAATCCCCATCCTCCATATATCCAAACAACAAAGCATAATATTTCGCCCACTAAGCCAATCACTTTATTGACTCCTAGCCGCAGACCTCCTCATTCTAACCTGAATCGGAGGACAACCAGTAAGCTACCCTTTTACCATCATTGGACAAGTAGCATCCGTACTATACTTCACAACAATCCTAATCCTAATACCAACTATCTCCCTAATTGAAAACAAAATACTCAAATGGGCCT")

mutation <- ("ATGACCCCAATACGCAAAACTAACCCCCTAATAAAATTAATTAACCACTCATTCATCGACCTCCCCACCCCATCCAACATCTCCGCATGATGAAACTTCGGCTCACTCCTTGGCGCCTGCCTGATCCTCCAAATCACCACAGGACTATTCCTAGCCATGCACTACTCACCAGACGCCTCAACCGCCTTTTCATCAATCGCCCACATCACTCGAGACGTAAATTATGGCTGAATCATCCGCTACCTTCACGCCAATGGCGCCTCAATATTCTTTATCTGCCTCTTCCTACACATCGGGCGAGGCCTATATTACGGATCATTTCTCTACTCAGAAACCTGAAACATCGGCATTATCCTCCTGCTTGCAACTATAGCAACAGCCTTCATAGGCTATGTCCTCCCGTGAGGCCAAATATCATTCTGAGGGGCCACAGTAATTACAAACTTACTATCCGCCATCCCATACATTGGGACAGACCTAGTTCAATGAATCTGAAGAGGCTACTCAGTAGACAGTCCCACCCTCACACGATTCTTTACCTTTCACTTCATCTTGCCCTTCATTATTGCAGCCCTAGCAACACTCCACCTCCTATTCTTGCACGAAACGGGATCAAACAACCCCCTAGGAATCACCTCCCATTCCGATAAAATCACCTTCCACCCTTACTACACAATCAAAGACGCCCTCGGCTTACTTCTCTTCCTTCTCTCCTTAATGACATTAACACTATTCTCACCAGACCTCCTAAGCGACCCAGACAATTATACCCTAGCCAACCCCTTAAACACCCCTCCCCACATCAAGCCCGAATGATATTTCCTATTCGCCTACACAATTCTCCGATCCGTCCCTAACAAACTAGGAGACGTCCTTGCCCTATTACTATCCATCCTCATCCTAGCAATAATCCCCATCCTCCATATATCCAAACAACAAAGCATAATATTTCGCCCACTAAGCCAATCACTTTATTGACTCCTAGCCGCAGACCTCCTCATTCTAACCTGAATCGGAGGACAACCAGTAAGCTACCCTTTTACCATCATTGGACAAGTAGCATCCGTACTATACTTCACAACAATCCTAATCCTAATACCAACTATCTCCCTAATTGAAAACAAAATACTCAAATGGGCCT")

library(diffobj)

diffObj(reference,mutation)

which.min(strsplit(reference, "")[[1]] == strsplit(mutation, "")[[1]])


difference <- function(a, b, exclude = c("-", "?"), ignore.case = TRUE, show.excluded = FALSE)
{
if(nchar(a)!=nchar(b)) stop("Lengths of input strings differ. Please check your input.")
if(ignore.case)
{
a <- toupper(a)
b <- toupper(b)
}
split_seqs <- strsplit(c(a, b), split = "")
only.diff <- (split_seqs[[1]] != split_seqs[[2]])
only.diff[
(split_seqs[[1]] %in% exclude) |
(split_seqs[[2]] %in% exclude)
] <- NA
diff.info<-data.frame(which(is.na(only.diff)|only.diff),
split_seqs[[1]][only.diff],split_seqs[[2]][only.diff])
names(diff.info)<-c("position","poly.seq.a","poly.seq.b")
if(!show.excluded) diff.info<-na.omit(diff.info)
diff.info
}

difference(reference,mutation)



reference1 <- substring(reference,
          seq(1, nchar(reference), 3),
          seq(n, nchar(reference),3))


ref <- as.data.frame(reference1)


mutation1 <- substring(mutation,
          seq(1, nchar(mutation), 3),
          seq(n, nchar(mutation),3))
mut <- as.data.frame(mutation1)

496/3 #165.3 (166 codon)
751/3 #250.3 (251 codon)
869/3 #289.6 (290 codon)

ref[166,] #Gly
mut[166,] #Arg

ref[251,] #Gly
mut[251,] #Ser

ref[290,] #Gly
mut[290,] #Asp


#############

library(readr)
library(dplyr)
library(stringi)

urlfile="https://raw.githubusercontent.com/NGSchoolEU/ngs22_registration_form/1cc647a3733e2c8a21b47aa497b4ca8c42457aa8/data/single-cell-studies.tsv"

scdata <- read.csv(url(urlfile), sep="\t", header=T)

scdata$Multiple<-sapply(stri_extract_all_regex(scdata$Organism,"[A-Z]+"),function(x) length(unique(x)))

scdata1 <- scdata[!(scdata$Multiple==1),]
nrow(scdata1)

scdata2 <- scdata[!(scdata$Multiple>1),]
nrow(scdata2)

y1 <- table(scdata2$Organism)
y2 <- y1[order(y1,decreasing = TRUE)]
head(y2)

colnames(scdata2)
cells <- as.data.frame(scdata2[,(8:9)])
cells$Reported.cells.total <- as.numeric(gsub(",","",cells$Reported.cells.total))
cells1 <- cells[order(cells$Reported.cells.total,decreasing=T),]
head(cells1)


colnames(scdata)

##tail(names(sort(table(organism2$Organism))), 1)


##test <- organism[order(organism$Multiple_Organisms,decreasing=T),]

0.109


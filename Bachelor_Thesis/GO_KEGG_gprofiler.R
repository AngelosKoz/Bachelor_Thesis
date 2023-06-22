setwd('/home/aggelos/Desktop/Rotation/Files')


##library("PloGO2")

library(org.Hs.eg.db)
library(parallel)
library(ggplot2)
library(gprofiler2)
library(tidyr)
library(plotly)
library(dplyr)

sessionInfo()

##protcounts <- read.delim("protein", sep="\t", header=T)
##protcounts <- as.data.frame(t(protcounts))
##protcounts <- protcounts[-c(1:3),]

##patients <- read.delim("patorderstages", sep="\t", header=T)
##patNOf0f1 <- as.data.frame(patients[cc(62:95),-c(2:5)])
##names(patNOf0f1) <- "patient"
##patNOf0f1



convnafl <- read.delim("controlnaflresults0,005adj", sep="\t", header=T)
convf0 <- read.delim("controlf0f1results0,005adj", sep="\t", header=T)
convf2 <- read.delim("controlf2results0,005adj", sep="\t", header=T)
convf3 <- read.delim("controlf3results0,005adj", sep="\t", header=T)
convf4 <- read.delim("controlf4results0,005adj", sep="\t", header=T)
naflvf2 <- read.delim("naflf2results0,005adj", sep="\t", header=T)
naflvf3 <- read.delim("naflf3results0,005adj", sep="\t", header=T)
naflvf4 <- read.delim("naflf4results0,005adj", sep="\t", header=T)
f0vf2 <- read.delim("f0f1vsf2results0,005adj", sep="\t", header=T)
f0vf3 <- read.delim("f0f1vsf3results0,005adj", sep="\t", header=T)
f0vf4 <- read.delim("f0f1vsf4results0,005adj", sep="\t", header=T)

convnaflorg <- read.delim("controlnaflresultsoriginal", sep="\t", header=T)
convf0org <- read.delim("controlf0f1resultsoriginal", sep="\t", header=T)
convf2org <- read.delim("controlf2resultsoriginal", sep="\t", header=T)
convf3org <- read.delim("controlf3resultsoriginal", sep="\t", header=T)
convf4org <- read.delim("controlf4resultsoriginal", sep="\t", header=T)
naflvf2org <- read.delim("naflf2resultsoriginal", sep="\t", header=T)
naflvf3org <- read.delim("naflf3resultsoriginal", sep="\t", header=T)
naflvf4org <- read.delim("naflf4resultsoriginal", sep="\t", header=T)
f0vf2org <- read.delim("f0f1vsf2resultsoriginal", sep="\t", header=T)
f0vf3org <- read.delim("f0f1vsf3resultsoriginal", sep="\t", header=T)
f0vf4org <- read.delim("f0f1vsf4resultsoriginal", sep="\t", header=T)


head(convnaflorg)
nrow(convnaflorg)

convnaflns1 <- convnaflorg[!is.na(convnaflorg$padj),]
convnaflns2 <- subset(convnaflns1, padj > 0.01, drop=FALSE) ##Confirm that not removing NA gives same result


convnaflns <- subset(convnaflorg, padj > 0.01, drop=FALSE) #11398
convf0ns <- subset(convf0org, padj > 0.01, drop=FALSE) #11633
convf2ns <- subset(convf2org, padj > 0.01, drop=FALSE) #11133
convf3ns <- subset(convf3org, padj > 0.01, drop=FALSE) #11760
convf4ns <- subset(convf4org, padj > 0.01, drop=FALSE) #11675
naflvf2ns <- subset(naflvf2org, padj > 0.01, drop=FALSE) #16265
naflvf3ns <- subset(naflvf3org, padj > 0.01, drop=FALSE) #14622
naflvf4ns <- subset(naflvf4org, padj > 0.01, drop=FALSE) #14710
f0vf2ns <- subset(f0vf2org, padj > 0.01, drop=FALSE) #15914
f0vf3ns <- subset(f0vf3org, padj > 0.01, drop=FALSE) #15416
f0vf4ns <- subset(f0vf4org, padj > 0.01, drop=FALSE) #14313

nrow(f0vf4ns)




##Control vs NAFL stage


head(convnafl)
nrow(convnafl) #4589

convnafl_UP <- subset(convnafl, log2FoldChange > 0.5849625) #1620
convnafl_UP <- as.data.frame(convnafl_UP[,-c(2:12)])
nrow(convnafl_UP)
names(convnafl_UP) <- "gene"

convnafl_UP_NEUT <- subset(convnafl, log2FoldChange >= 0) #2445
convnafl_UP_NEUT <- as.data.frame(convnafl_UP_NEUT[,-c(2:12)])
nrow(convnafl_UP_NEUT)
names(convnafl_UP_NEUT) <- "gene"

convnafl_DOWN <- subset(convnafl, log2FoldChange < -0.5849625) #1554
convnafl_DOWN <- as.data.frame(convnafl_DOWN[,-c(2:12)])
nrow(convnafl_DOWN)
names(convnafl_DOWN) <- "gene"

convnafl_DOWN_NEUT <- subset(convnafl, log2FoldChange <= 0) #2144
convnafl_DOWN_NEUT <- as.data.frame(convnafl_DOWN_NEUT[,-c(2:12)])
nrow(convnafl_DOWN_NEUT)
names(convnafl_DOWN_NEUT) <- "gene"




###Control vs NASH F0-F1

nrow(convf0) #4238

convf0_UP <- subset(convf0, log2FoldChange > 0.5849625) #1292
convf0_UP <- as.data.frame(convf0_UP[,-c(2:12)])
nrow(convf0_UP)
names(convf0_UP) <- "gene"

convf0_UP_NEUT <- subset(convf0, log2FoldChange >= 0) #2077
convf0_UP_NEUT <- as.data.frame(convf0_UP_NEUT[,-c(2:12)])
nrow(convf0_UP_NEUT)
names(convf0_UP_NEUT) <- "gene"

convf0_DOWN <- subset(convf0, log2FoldChange < -0.5849625) #1547
convf0_DOWN <- as.data.frame(convf0_DOWN[,-c(2:12)])
nrow(convf0_DOWN)
names(convf0_DOWN) <- "gene"

convf0_DOWN_NEUT <- subset(convf0, log2FoldChange <= 0) #2161
convf0_DOWN_NEUT <- as.data.frame(convf0_DOWN_NEUT[,-c(2:12)])
nrow(convf0_DOWN_NEUT)
names(convf0_DOWN_NEUT) <- "gene"


###Control vs NASH F2

nrow(convf2) #4853

convf2_UP <- subset(convf2, log2FoldChange > 0.5849625) #1815
convf2_UP <- as.data.frame(convf2_UP[,-c(2:12)])
 nrow(convf2_UP)
names(convf2_UP) <- "gene"

convf2_UP_NEUT <- subset(convf2, log2FoldChange >= 0) #2602
convf2_UP_NEUT <- as.data.frame(convf2_UP_NEUT[,-c(2:12)])
nrow(convf2_UP_NEUT)
names(convf2_UP_NEUT) <- "gene"

convf2_DOWN <- subset(convf2, log2FoldChange < -0.5849625) #1701
convf2_DOWN <- as.data.frame(convf2_DOWN[,-c(2:12)])
nrow(convf2_DOWN)
names(convf2_DOWN) <- "gene"

convf2_DOWN_NEUT <- subset(convf2, log2FoldChange <= 0) #2251
convf2_DOWN_NEUT <- as.data.frame(convf2_DOWN_NEUT[,-c(2:12)])
nrow(convf2_DOWN_NEUT)
names(convf2_DOWN_NEUT) <- "gene"


###Control vs NASH F3

nrow(convf3) #4275

convf3_UP <- subset(convf3, log2FoldChange > 0.5849625) #1488
convf3_UP <- as.data.frame(convf3_UP[,-c(2:12)])
nrow(convf3_UP)
names(convf3_UP) <- "gene"

convf3_UP_NEUT <- subset(convf3, log2FoldChange >= 0) #2238
convf3_UP_NEUT <- as.data.frame(convf3_UP_NEUT[,-c(2:12)])
nrow(convf3_UP_NEUT)
names(convf3_UP_NEUT) <- "gene"

convf3_DOWN <- subset(convf3, log2FoldChange < -0.5849625) #1368
convf3_DOWN <- as.data.frame(convf3_DOWN[,-c(2:12)])
nrow(convf3_DOWN)
names(convf3_DOWN) <- "gene"

convf3_DOWN_NEUT <- subset(convf3, log2FoldChange <= 0) #2037
convf3_DOWN_NEUT <- as.data.frame(convf3_DOWN_NEUT[,-c(2:12)])
nrow(convf3_DOWN_NEUT)
names(convf3_DOWN_NEUT) <- "gene"

###Control vs NASHF4

nrow(convf4) #3661

convf4_UP <- subset(convf4, log2FoldChange > 0.5849625) #1599
convf4_UP <- as.data.frame(convf4_UP[,-c(2:12)])
nrow(convf4_UP)
names(convf4_UP) <- "gene"

convf4_UP_NEUT <- subset(convf4, log2FoldChange >= 0) #1926
convf4_UP_NEUT <- as.data.frame(convf4_UP_NEUT[,-c(2:12)])
nrow(convf4_UP_NEUT)
names(convf4_UP_NEUT) <- "gene"

convf4_DOWN <- subset(convf4, log2FoldChange < -0.5849625) #1440
convf4_DOWN <- as.data.frame(convf4_DOWN[,-c(2:12)])
nrow(convf4_DOWN)
names(convf4_DOWN) <- "gene"

convf4_DOWN_NEUT <- subset(convf4, log2FoldChange <= 0) #1735
convf4_DOWN_NEUT <- as.data.frame(convf4_DOWN_NEUT[,-c(2:12)])
nrow(convf4_DOWN_NEUT)
names(convf4_DOWN_NEUT) <- "gene"

###NAFL vs NASH F2

nrow(naflvf2)#58

naflvf2_UP <- subset(naflvf2, log2FoldChange > 0.5849625) #43
naflvf2_UP <- as.data.frame(naflvf2_UP[,-c(2:12)])
nrow(naflvf2_UP)
names(naflvf2_UP) <- "gene"

naflvf2_UP_NEUT <- subset(naflvf2, log2FoldChange >= 0) #48
naflvf2_UP_NEUT <- as.data.frame(naflvf2_UP_NEUT[,-c(2:12)])
nrow(naflvf2_UP_NEUT)
names(naflvf2_UP_NEUT) <- "gene"

naflvf2_DOWN <- subset(naflvf2, log2FoldChange < -0.5849625) #2
naflvf2_DOWN <- as.data.frame(naflvf2_DOWN[,-c(2:12)])
nrow(naflvf2_DOWN)
names(naflvf2_DOWN) <- "gene"

naflvf2_DOWN_NEUT <- subset(naflvf2, log2FoldChange <= 0) #10
naflvf2_DOWN_NEUT <- as.data.frame(naflvf2_DOWN_NEUT[,-c(2:12)])
nrow(naflvf2_DOWN_NEUT)
names(naflvf2_DOWN_NEUT) <- "gene"



### NAFL vs NASH F3

nrow(naflvf3)#1389

naflvf3_UP <- subset(naflvf3, log2FoldChange > 0.5849625) #668
naflvf3_UP <- as.data.frame(naflvf3_UP[,-c(2:12)])
nrow(naflvf3_UP)
names(naflvf3_UP) <- "gene"

naflvf3_UP_NEUT <- subset(naflvf3, log2FoldChange >= 0) #929
naflvf3_UP_NEUT <- as.data.frame(naflvf3_UP_NEUT[,-c(2:12)])
nrow(naflvf3_UP_NEUT)
names(naflvf3_UP_NEUT) <- "gene"

naflvf3_DOWN <- subset(naflvf3, log2FoldChange < -0.5849625) #70
naflvf3_DOWN <- as.data.frame(naflvf3_DOWN[,-c(2:12)])
nrow(naflvf3_DOWN)
names(naflvf3_DOWN) <- "gene"

naflvf3_DOWN_NEUT <- subset(naflvf3, log2FoldChange <= 0) #460
naflvf3_DOWN_NEUT <- as.data.frame(naflvf3_DOWN_NEUT[,-c(2:12)])
nrow(naflvf3_DOWN_NEUT)
names(naflvf3_DOWN_NEUT) <- "gene"

### NAFL vs NASH F4

nrow(naflvf4) #1551

naflvf4_UP <- subset(naflvf4, log2FoldChange > 0.5849625) #980
naflvf4_UP <- as.data.frame(naflvf4_UP[,-c(2:12)])
nrow(naflvf4_UP)
names(naflvf4_UP) <- "gene"

naflvf4_UP_NEUT <- subset(naflvf4, log2FoldChange >= 0) #1103
naflvf4_UP_NEUT <- as.data.frame(naflvf4_UP_NEUT[,-c(2:12)])
nrow(naflvf4_UP_NEUT)
names(naflvf4_UP_NEUT) <- "gene"

naflvf4_DOWN <- subset(naflvf4, log2FoldChange < -0.5849625) #159
naflvf4_DOWN <- as.data.frame(naflvf4_DOWN[,-c(2:12)])
nrow(naflvf4_DOWN)
names(naflvf4_DOWN) <- "gene"

naflvf4_DOWN_NEUT <- subset(naflvf4, log2FoldChange <= 0) #448
naflvf4_DOWN_NEUT <- as.data.frame(naflvf4_DOWN_NEUT[,-c(2:12)])
nrow(naflvf4_DOWN_NEUT)
names(naflvf4_DOWN_NEUT) <- "gene"

### NASH F0-F1 vs NASH F2


f0vf2

nrow(f0vf2) #11

f0vf2_UP <- subset(f0vf2, log2FoldChange > 0.5849625) #11
f0vf2_UP <- as.data.frame(f0vf2_UP[,-c(2:12)])
nrow(f0vf2_UP)
names(f0vf2_UP) <- "gene"

f0vf2_UP_NEUT <- subset(f0vf2, log2FoldChange >= 0) #11
f0vf2_UP_NEUT <- as.data.frame(f0vf2_UP_NEUT[,-c(2:12)])
nrow(f0vf2_UP_NEUT)
names(f0vf2_UP_NEUT) <- "gene"

f0vf2_DOWN <- subset(f0vf2, log2FoldChange < -0.5849625) #0
f0vf2_DOWN <- as.data.frame(f0vf2_DOWN[,-c(2:12)])
nrow(f0vf2_DOWN)
names(f0vf2_DOWN) <- "gene"

f0vf2_DOWN_NEUT <- subset(f0vf2, log2FoldChange <= 0) #0
f0vf2_DOWN_NEUT <- as.data.frame(f0vf2_DOWN_NEUT[,-c(2:12)])
nrow(f0vf2_DOWN_NEUT)
names(f0vf2_DOWN_NEUT) <- "gene"

### NASH F0-F1 vs NASH F3

nrow(f0vf3) #756

f0vf3_UP <- subset(f0vf3, log2FoldChange > 0.5849625) #556
f0vf3_UP <- as.data.frame(f0vf3_UP[,-c(2:12)])
nrow(f0vf3_UP)
names(f0vf3_UP) <- "gene"

f0vf3_UP_NEUT <- subset(f0vf3, log2FoldChange >= 0) #668
f0vf3_UP_NEUT <- as.data.frame(f0vf3_UP_NEUT[,-c(2:12)])
nrow(f0vf3_UP_NEUT)
names(f0vf3_UP_NEUT) <- "gene"

f0vf3_DOWN <- subset(f0vf3, log2FoldChange < -0.5849625) #25
f0vf3_DOWN <- as.data.frame(f0vf3_DOWN[,-c(2:12)])
nrow(f0vf3_DOWN)
names(f0vf3_DOWN) <- "gene"

f0vf3_DOWN_NEUT <- subset(f0vf3, log2FoldChange <= 0) #88
f0vf3_DOWN_NEUT <- as.data.frame(f0vf3_DOWN_NEUT[,-c(2:12)])
nrow(f0vf3_DOWN_NEUT)
names(f0vf3_DOWN_NEUT) <- "gene"


### NASH F0-F1 vs NASH F4

nrow(f0vf4) #1769

f0vf4_UP <- subset(f0vf4, log2FoldChange > 0.5849625) #1191
f0vf4_UP <- as.data.frame(f0vf4_UP[,-c(2:12)])
nrow(f0vf4_UP)
names(f0vf4_UP) <- "gene"

f0vf4_UP_NEUT <- subset(f0vf4, log2FoldChange >= 0) #1320
f0vf4_UP_NEUT <- as.data.frame(f0vf4_UP_NEUT[,-c(2:12)])
nrow(f0vf4_UP_NEUT)
names(f0vf4_UP_NEUT) <- "gene"

f0vf4_DOWN <- subset(f0vf4, log2FoldChange < -0.5849625) #173
f0vf4_DOWN <- as.data.frame(f0vf4_DOWN[,-c(2:12)])
nrow(f0vf4_DOWN)
names(f0vf4_DOWN) <- "gene"

f0vf4_DOWN_NEUT <- subset(f0vf4, log2FoldChange <= 0) #449
f0vf4_DOWN_NEUT <- as.data.frame(f0vf4_DOWN_NEUT[,-c(2:12)])
nrow(f0vf4_DOWN_NEUT)
names(f0vf4_DOWN_NEUT) <- "gene"



#geneid <- as.data.frame(convnafl[,1]) ##CHANGE when you input new set of genes
#names(geneid) <- "gene"
#geneid

geneid1up <- convnafl_UP
geneid2up <- convf0_UP
geneid3up <- convf2_UP
geneid4up <-convf3_UP
geneid5up <-convf4_UP
geneid6up <-naflvf2_UP
geneid7up <-naflvf3_UP
geneid8up <-naflvf4_UP
geneid9up <-f0vf2_UP
geneid10up <-f0vf3_UP
geneid11up <-f0vf4_UP

#geneid1upns <- convnafl_UP_NS
#geneid2upns <- convf0_UP_NS
#geneid3upns <- convf2_UP_NS
#geneid4upns <-convf3_UP_NS
#geneid5upns <-convf4_UP_NS
#geneid6upns <-naflvf2_UP_NS
#geneid7upns <-naflvf3_UP_NS
#geneid8upns <-naflvf4_UP_NS
#geneid9upns <-f0vf2_UP_NS
#geneid10upns <-f0vf3_UP_NS
#geneid11upns <-f0vf4_UP_NS

geneid1upneut <- convnafl_UP_NEUT
geneid2upneut <- convf0_UP_NEUT
geneid3upneut <- convf2_UP_NEUT
geneid4upneut <-convf3_UP_NEUT
geneid5upneut <-convf4_UP_NEUT
geneid6upneut <-naflvf2_UP_NEUT
geneid7upneut <-naflvf3_UP_NEUT
geneid8upneut <-naflvf4_UP_NEUT
geneid9upneut <-f0vf2_UP_NEUT
geneid10upneut <-f0vf3_UP_NEUT
geneid11upneut <-f0vf4_UP_NEUT

geneid1down <- convnafl_DOWN
geneid2down <- convf0_DOWN
geneid3down <- convf2_DOWN
geneid4down <-convf3_DOWN
geneid5down <-convf4_DOWN
geneid6down <-naflvf2_DOWN
geneid7down <-naflvf3_DOWN
geneid8down <-naflvf4_DOWN
geneid10down <-f0vf3_DOWN
geneid11down <-f0vf4_DOWN

geneid1downneut <- convnafl_DOWN_NEUT
geneid2downneut <- convf0_DOWN_NEUT
geneid3downneut <- convf2_DOWN_NEUT
geneid4downneut <-convf3_DOWN_NEUT
geneid5downneut <-convf4_DOWN_NEUT
geneid6downneut <-naflvf2_DOWN_NEUT
geneid7downneut <-naflvf3_DOWN_NEUT
geneid8downneut <-naflvf4_DOWN_NEUT
geneid10downneut <-f0vf3_DOWN_NEUT
geneid11downneut <-f0vf4_DOWN_NEUT

#geneid1downns <- convnafl_DOWN_NS
#geneid2downns <- convf0_DOWN_NS
#geneid3downns <- convf2_DOWN_NS
#geneid4downns <- convf3_DOWN_NS
#geneid5downns <- convf4_DOWN_NS
#geneid6downns <- naflvf2_DOWN_NS
#geneid7downns <- naflvf3_DOWN_NS
#geneid8downns <- naflvf4_DOWN_NS
#geneid9downns <- f0vf2_DOWN_NS
#geneid10downns <- f0vf3_DOWN_NS
#geneid11downns <- f0vf4_DOWN_NS


                                        #Gene Ontology Analysis







## gProfiler GO/KEGG



GOresult <- gost(
geneid1up$gene,   ###Try test after gconvert//geneid1 originally
organism = "hsapiens",
ordered_query = FALSE, ##Order based on FC and change this to TRUE??
multi_query = FALSE,
significant = TRUE, 
exclude_iea = FALSE,
measure_underrepresentation = FALSE,
evcodes = FALSE, #Change to true
user_threshold = 0.05,
correction_method =  "gSCS", ## c("g_SCS", "bonferroni", "fdr", "false_discovery_rate","analytical") Different correction methods
domain_scope = "annotated", ## c("known", "custom", "custom_annotated")
custom_bg = NULL,
numeric_ns = "",
sources = c("GO:BP","GO:MF","GO:CC","KEGG"),
as_short_link = FALSE)

GOresult1 <- as.data.frame(GOresult$result)
GOresult1$minuslog10pval <- -log10(GOresult1$p_value)
names(GOresult1)[15] <- "-log10(pval)"
head(GOresult1)

write.table(GOresult1, "ConVnafl_GO_KEGG_UP_0.005", quote=F, sep="\t")

GOresult2 <- GOresult1[order(GOresult1$p_value, decreasing=F),]
head(GOresult2)

GOCC <- GOresult2[GOresult2$source == "GO:BP",]
GOMF <- GOresult2[GOresult2$source == "GO:MF",]
GOBP <- GOresult2[GOresult2$source == "GO:CC",]
KEGG <- GOresult2[GOresult2$source == "KEGG",]

rownames(GOBP) <- 1:nrow(GOBP)
head(GOBP)
rownames(GOMF) <- 1:nrow(GOMF)
head(GOMF)
rownames(GOCC) <- 1:nrow(GOCC)
head(GOCC)
rownames(KEGG) <- 1:nrow(KEGG)
head(KEGG)

                                        #Full analysis Plot


publish_gosttable(GOresult2, highlight_terms = GOresult2[c(1:20),],
                        use_colors = FALSE, 
                        show_columns = c("source", "term_name", "term_size", "intersection_size","-log10(pval)"),
                  filename = NULL)+
    ggtitle("Mytest")




ggsave("test.jpeg")


                                        #BP analysis Plot

plot2 <- publish_gosttable(GOBP, highlight_terms = GOBP[c(1:20),],
                        use_colors = FALSE, 
                        show_columns = c("source", "term_name", "term_size", "intersection_size","-log10(pval)"),
                        filename = NULL)
plot2

ggsave("test.jpeg")

                                        #MF analysis Plot

plot3 <- publish_gosttable(GOMF, highlight_terms = GOMF[c(1:20),],
                        use_colors = FALSE, 
                        show_columns = c("source", "term_name", "term_size", "intersection_size","-log10(pval)"),
                        filename = NULL)
plot3

ggsave("test.jpeg")


                                        #CC analysis Plot

plot4 <- publish_gosttable(GOCC, highlight_terms = GOCC[c(1:20),],
                        use_colors = FALSE, 
                        show_columns = c("source", "term_name", "term_size", "intersection_size","-log10(pval)"),
                        filename = NULL)
plot4

ggsave("test.jpeg")


                                        #KEGG analysis Plot

plot5 <- publish_gosttable(KEGG, highlight_terms = KEGG[c(1:20),],
                        use_colors = FALSE, 
                        show_columns = c("source", "term_name", "term_size", "intersection_size","-log10(pval)"),
                        filename = NULL)
plot5

ggsave("test.jpeg")


                                        #OTHER FUNCTIONS AND LEGENDS


##################
gconvert ##gene/protein/transcript identifier conversion
gorth ##Orthology Search across species
gsnpense ##mapping SNP rs-identifiers to chromosome positions, protein coding genes and variant effects
##################


###QUERY - the name of the input query which by default is the order of query with the prefix “query_.” This can be changed by using a named list input.
###TERM_SIZE - number of genes that are annotated to the term
###QUERY_SIZE - number of genes that were included in the query. This might be different from the size of the original list if:
#1)any genes were mapped to multiple Ensembl gene IDs
#2)any genes failed to be mapped to Ensembl gene IDs
#3)the parameter ordered_query = TRUE and the optimal cutoff for the term was found before the end of the query
#4)the domain_scope was set to “annotated” or “custom”
###INTERSECTION_SIZE - the number of genes in the input query that are annotated to the corresponding term
###PRECISION - the proportion of genes in the input list that are annotated to the function (defined as intersection_size/query_size)
###RECALL - the proportion of functionally annotated genes that the query recovers (defined as intersection_size/term_size)
###TERM_ID - unique term identifier (e.g GO:0005005)
###SOURCE - the abbreviation of the data source for the term (e.g. GO:BP)
###TERM_NAME - the short name of the function
###EFFECTIVE_DOMAIN_SIZE - the total number of genes “in the universe” used for the hypergeometric test
###SOURCE_ORDER - numeric order for the term within its data source (this is important for drawing the results)
###PARENTS - list of term IDs that are hierarchically directly above the term. For non-hierarchical data sources this points to an artificial root node.

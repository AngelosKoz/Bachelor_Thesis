setwd('/home/aggelos/Desktop/Rotation/Files')


##library("PloGO2")

library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)

library(parallel)

library(ggplot2)

sessionInfo()

##protcounts <- read.delim("protein", sep="\t", header=T)
##protcounts <- as.data.frame(t(protcounts))
##protcounts <- protcounts[-c(1:3),]

##patients <- read.delim("patorderstages", sep="\t", header=T)
##patNOf0f1 <- as.data.frame(patients[cc(62:95),-c(2:5)])
##names(patNOf0f1) <- "patient"
##patNOf0f1


###### Do the same with 0,005/0,001 + LNC

convnafl <- read.delim("controlnaflresults0,01adj", sep="\t", header=T)
convf0 <- read.delim("controlf0f1results0,01adj", sep="\t", header=T)
convf2 <- read.delim("controlf2results0,01adj", sep="\t", header=T)
convf3 <- read.delim("controlf3results0,01adj", sep="\t", header=T)
convf4 <- read.delim("controlf4results0,01adj", sep="\t", header=T)
naflvf2 <- read.delim("naflf2results0,01adj", sep="\t", header=T)
naflvf3 <- read.delim("naflf3results0,01adj", sep="\t", header=T)
naflvf4 <- read.delim("naflf4results0,01adj", sep="\t", header=T)
f0vf2 <- read.delim("f0f1vsf2results0,01adj", sep="\t", header=T)
f0vf3 <- read.delim("f0f1vsf3results0,01adj", sep="\t", header=T)
f0vf4 <- read.delim("f0f1vsf4results0,01adj", sep="\t", header=T)

##Control vs NAFL stage


head(convnafl)
nrow(convnafl) #5135

convnafl_UP <- subset(convnafl, log2FoldChange > 0.5849625) #1728
convnafl_UP <- as.data.frame(convnafl_UP[,-c(2:12)])
nrow(convnafl_UP)
names(convnafl_UP) <- "gene"

convnafl_UP_NEUT <- subset(convnafl, log2FoldChange >= 0) #2679
convnafl_UP_NEUT <- as.data.frame(convnafl_UP_NEUT[,-c(2:12)])
nrow(convnafl_UP_NEUT)
names(convnafl_UP_NEUT) <- "gene"

convnafl_DOWN <- subset(convnafl, log2FoldChange < -0.5849625) #1716
convnafl_DOWN <- as.data.frame(convnafl_DOWN[,-c(2:12)])
nrow(convnafl_DOWN)
names(convnafl_DOWN) <- "gene"

convnafl_DOWN_NEUT <- subset(convnafl, log2FoldChange <= 0) #2456
convnafl_DOWN_NEUT <- as.data.frame(convnafl_DOWN_NEUT[,-c(2:12)])
nrow(convnafl_DOWN_NEUT)
names(convnafl_DOWN_NEUT) <- "gene"


###Control vs NASH F0-F1

nrow(convf0) #4771

convf0_UP <- subset(convf0, log2FoldChange > 0.5849625) #1385
convf0_UP <- as.data.frame(convf0_UP[,-c(2:12)])
nrow(convf0_UP)
names(convf0_UP) <- "gene"

convf0_UP_NEUT <- subset(convf0, log2FoldChange >= 0) #2322
convf0_UP_NEUT <- as.data.frame(convf0_UP_NEUT[,-c(2:12)])
nrow(convf0_UP_NEUT)
names(convf0_UP_NEUT) <- "gene"

convf0_DOWN <- subset(convf0, log2FoldChange < -0.5849625) #1665
convf0_DOWN <- as.data.frame(convf0_DOWN[,-c(2:12)])
nrow(convf0_DOWN)
names(convf0_DOWN) <- "gene"

convf0_DOWN_NEUT <- subset(convf0, log2FoldChange <= 0) #2449
convf0_DOWN_NEUT <- as.data.frame(convf0_DOWN_NEUT[,-c(2:12)])
nrow(convf0_DOWN_NEUT)
names(convf0_DOWN_NEUT) <- "gene"


###Control vs NASH F2

nrow(convf2) #5417

convf2_UP <- subset(convf2, log2FoldChange > 0.5849625) #1948
convf2_UP <- as.data.frame(convf2_UP[,-c(2:12)])
nrow(convf2_UP)
names(convf2_UP) <- "gene"

convf2_UP_NEUT <- subset(convf2, log2FoldChange >= 0) #2859
convf2_UP_NEUT <- as.data.frame(convf2_UP_NEUT[,-c(2:12)])
nrow(convf2_UP_NEUT)
names(convf2_UP_NEUT) <- "gene"

convf2_DOWN <- subset(convf2, log2FoldChange < -0.5849625) #1819
convf2_DOWN <- as.data.frame(convf2_DOWN[,-c(2:12)])
nrow(convf2_DOWN)
names(convf2_DOWN) <- "gene"

convf2_DOWN_NEUT <- subset(convf2, log2FoldChange <= 0) #2558
convf2_DOWN_NEUT <- as.data.frame(convf2_DOWN_NEUT[,-c(2:12)])
nrow(convf2_DOWN_NEUT)
names(convf2_DOWN_NEUT) <- "gene"


###Control vs NASH F3

nrow(convf3) #4854

convf3_UP <- subset(convf3, log2FoldChange > 0.5849625) #1645
convf3_UP <- as.data.frame(convf3_UP[,-c(2:12)])
nrow(convf3_UP)
names(convf3_UP) <- "gene"

convf3_UP_NEUT <- subset(convf3, log2FoldChange >= 0) #2546
convf3_UP_NEUT <- as.data.frame(convf3_UP_NEUT[,-c(2:12)])
nrow(convf3_UP_NEUT)
names(convf3_UP_NEUT) <- "gene"

convf3_DOWN <- subset(convf3, log2FoldChange < -0.5849625) #1474
convf3_DOWN <- as.data.frame(convf3_DOWN[,-c(2:12)])
nrow(convf3_DOWN)
names(convf3_DOWN) <- "gene"

convf3_DOWN_NEUT <- subset(convf3, log2FoldChange <= 0) #2308
convf3_DOWN_NEUT <- as.data.frame(convf3_DOWN_NEUT[,-c(2:12)])
nrow(convf3_DOWN_NEUT)
names(convf3_DOWN_NEUT) <- "gene"

###Control vs NASHF4

nrow(convf4) #4259

convf4_UP <- subset(convf4, log2FoldChange > 0.5849625) #1813
convf4_UP <- as.data.frame(convf4_UP[,-c(2:12)])
nrow(convf4_UP)
names(convf4_UP) <- "gene"

convf4_UP_NEUT <- subset(convf4, log2FoldChange >= 0) #2243
convf4_UP_NEUT <- as.data.frame(convf4_UP_NEUT[,-c(2:12)])
nrow(convf4_UP_NEUT)
names(convf4_UP_NEUT) <- "gene"

convf4_DOWN <- subset(convf4, log2FoldChange < -0.5849625) #1589
convf4_DOWN <- as.data.frame(convf4_DOWN[,-c(2:12)])
nrow(convf4_DOWN)
names(convf4_DOWN) <- "gene"

convf4_DOWN_NEUT <- subset(convf4, log2FoldChange <= 0) #2016
convf4_DOWN_NEUT <- as.data.frame(convf4_DOWN_NEUT[,-c(2:12)])
nrow(convf4_DOWN_NEUT)
names(convf4_DOWN_NEUT) <- "gene"

###NAFL vs NASH F2

nrow(naflvf2)95

naflvf2_UP <- subset(naflvf2, log2FoldChange > 0.5849625) #59
naflvf2_UP <- as.data.frame(naflvf2_UP[,-c(2:12)])
nrow(naflvf2_UP)
names(naflvf2_UP) <- "gene"

naflvf2_UP_NEUT <- subset(naflvf2, log2FoldChange >= 0) #74
naflvf2_UP_NEUT <- as.data.frame(naflvf2_UP_NEUT[,-c(2:12)])
nrow(naflvf2_UP_NEUT)
names(naflvf2_UP_NEUT) <- "gene"

naflvf2_DOWN <- subset(naflvf2, log2FoldChange < -0.5849625) #4
naflvf2_DOWN <- as.data.frame(naflvf2_DOWN[,-c(2:12)])
nrow(naflvf2_DOWN)
names(naflvf2_DOWN) <- "gene"

naflvf2_DOWN_NEUT <- subset(naflvf2, log2FoldChange <= 0) #21
naflvf2_DOWN_NEUT <- as.data.frame(naflvf2_DOWN_NEUT[,-c(2:12)])
nrow(naflvf2_DOWN_NEUT)
names(naflvf2_DOWN_NEUT) <- "gene"



### NAFL vs NASH F3

nrow(naflvf3)#1791

naflvf3_UP <- subset(naflvf3, log2FoldChange > 0.5849625) #768
naflvf3_UP <- as.data.frame(naflvf3_UP[,-c(2:12)])
nrow(naflvf3_UP)
names(naflvf3_UP) <- "gene"

naflvf3_UP_NEUT <- subset(naflvf3, log2FoldChange >= 0) #1149
naflvf3_UP_NEUT <- as.data.frame(naflvf3_UP_NEUT[,-c(2:12)])
nrow(naflvf3_UP_NEUT)
names(naflvf3_UP_NEUT) <- "gene"

naflvf3_DOWN <- subset(naflvf3, log2FoldChange < -0.5849625) #79
naflvf3_DOWN <- as.data.frame(naflvf3_DOWN[,-c(2:12)])
nrow(naflvf3_DOWN)
names(naflvf3_DOWN) <- "gene"

naflvf3_DOWN_NEUT <- subset(naflvf3, log2FoldChange <= 0) #642
naflvf3_DOWN_NEUT <- as.data.frame(naflvf3_DOWN_NEUT[,-c(2:12)])
nrow(naflvf3_DOWN_NEUT)
names(naflvf3_DOWN_NEUT) <- "gene"

### NAFL vs NASH F4

nrow(naflvf4) #1886

naflvf4_UP <- subset(naflvf4, log2FoldChange > 0.5849625) #1118
naflvf4_UP <- as.data.frame(naflvf4_UP[,-c(2:12)])
nrow(naflvf4_UP)
names(naflvf4_UP) <- "gene"

naflvf4_UP_NEUT <- subset(naflvf4, log2FoldChange >= 0) #1304
naflvf4_UP_NEUT <- as.data.frame(naflvf4_UP_NEUT[,-c(2:12)])
nrow(naflvf4_UP_NEUT)
names(naflvf4_UP_NEUT) <- "gene"

naflvf4_DOWN <- subset(naflvf4, log2FoldChange < -0.5849625) #192
naflvf4_DOWN <- as.data.frame(naflvf4_DOWN[,-c(2:12)])
nrow(naflvf4_DOWN)
names(naflvf4_DOWN) <- "gene"

naflvf4_DOWN_NEUT <- subset(naflvf4, log2FoldChange <= 0) #582
naflvf4_DOWN_NEUT <- as.data.frame(naflvf4_DOWN_NEUT[,-c(2:12)])
nrow(naflvf4_DOWN_NEUT)
names(naflvf4_DOWN_NEUT) <- "gene"

### NASH F0-F1 vs NASH F2

nrow(f0vf2) #25

f0vf2_UP <- subset(f0vf2, log2FoldChange > 0.5849625) #21
f0vf2_UP <- as.data.frame(f0vf2_UP[,-c(2:12)])
nrow(f0vf2_UP)
names(f0vf2_UP) <- "gene"

f0vf2_UP_NEUT <- subset(f0vf2, log2FoldChange >= 0) #23
f0vf2_UP_NEUT <- as.data.frame(f0vf2_UP_NEUT[,-c(2:12)])
nrow(f0vf2_UP_NEUT)
names(f0vf2_UP_NEUT) <- "gene"

f0vf2_DOWN <- subset(f0vf2, log2FoldChange < -0.5849625) #0
f0vf2_DOWN <- as.data.frame(f0vf2_DOWN[,-c(2:12)])
nrow(f0vf2_DOWN)
names(f0vf2_DOWN) <- "gene"

f0vf2_DOWN_NEUT <- subset(f0vf2, log2FoldChange <= 0) #2
f0vf2_DOWN_NEUT <- as.data.frame(f0vf2_DOWN_NEUT[,-c(2:12)])
nrow(f0vf2_DOWN_NEUT)
names(f0vf2_DOWN_NEUT) <- "gene"

### NASH F0-F1 vs NASH F3

nrow(f0vf3) #924

f0vf3_UP <- subset(f0vf3, log2FoldChange > 0.5849625) #626
f0vf3_UP <- as.data.frame(f0vf3_UP[,-c(2:12)])
nrow(f0vf3_UP)
names(f0vf3_UP) <- "gene"

f0vf3_UP_NEUT <- subset(f0vf3, log2FoldChange >= 0) #791
f0vf3_UP_NEUT <- as.data.frame(f0vf3_UP_NEUT[,-c(2:12)])
nrow(f0vf3_UP_NEUT)
names(f0vf3_UP_NEUT) <- "gene"

f0vf3_DOWN <- subset(f0vf3, log2FoldChange < -0.5849625) #32
f0vf3_DOWN <- as.data.frame(f0vf3_DOWN[,-c(2:12)])
nrow(f0vf3_DOWN)
names(f0vf3_DOWN) <- "gene"

f0vf3_DOWN_NEUT <- subset(f0vf3, log2FoldChange <= 0) #133
f0vf3_DOWN_NEUT <- as.data.frame(f0vf3_DOWN_NEUT[,-c(2:12)])
nrow(f0vf3_DOWN_NEUT)
names(f0vf3_DOWN_NEUT) <- "gene"


### NASH F0-F1 vs NASH F4

nrow(f0vf4) #2167

f0vf4_UP <- subset(f0vf4, log2FoldChange > 0.5849625) #1354
f0vf4_UP <- as.data.frame(f0vf4_UP[,-c(2:12)])
nrow(f0vf4_UP)
names(f0vf4_UP) <- "gene"

f0vf4_UP_NEUT <- subset(f0vf4, log2FoldChange >= 0) #1563
f0vf4_UP_NEUT <- as.data.frame(f0vf4_UP_NEUT[,-c(2:12)])
nrow(f0vf4_UP_NEUT)
names(f0vf4_UP_NEUT) <- "gene"

f0vf4_DOWN <- subset(f0vf4, log2FoldChange < -0.5849625) #204
f0vf4_DOWN <- as.data.frame(f0vf4_DOWN[,-c(2:12)])
nrow(f0vf4_DOWN)
names(f0vf4_DOWN) <- "gene"

f0vf4_DOWN_NEUT <- subset(f0vf4, log2FoldChange <= 0) #604
f0vf4_DOWN_NEUT <- as.data.frame(f0vf4_DOWN_NEUT[,-c(2:12)])
nrow(f0vf4_DOWN_NEUT)
names(f0vf4_DOWN_NEUT) <- "gene"

##Unecessary part for now
#test <- read.delim("MLproteindf", sep="\t", header=T) 
#test1 <- test[-c(1:10,62:95), -c(18402:18405)] ##Remove patients not used 
#test2 <- as.data.frame(t(test1))
#test3 <- tibble::rownames_to_column(test2, "gene")
#test4 <- merge(geneid, test3, by="gene")


geneid <- as.data.frame(convnafl[,1]) ##CHANGE when you input new set of genes
names(geneid) <- "gene"
geneid

geneid1 <- convnafl_UP
geneid2 <- convf0_UP
geneid3 <- convf2_UP
geneid4 <-convf3_UP
geneid5 <-convf4_UP
geneid6 <-naflvf2_UP
geneid7 <-naflvf3_UP
geneid8 <-naflvf4_UP
geneid9 <-f0vf2_UP
geneid10 <-f0vf3_UP
geneid11 <-f0vf4_UP

geneid12 <-convnafl_UP_NEUT
geneid13 <-convf0_UP_NEUT
geneid14 <-convf2_UP_NEUT
geneid15 <-convf3_UP_NEUT
geneid16 <-convf4_UP_NEUT
geneid17 <-naflvf2_UP_NEUT
geneid18 <-naflvf3_UP_NEUT
geneid19 <-naflvf4_UP_NEUT
geneid20 <-f0vf2_UP_NEUT
geneid21<-f0vf3_UP_NEUT
geneid22 <-f0vf4_UP_NEUT

geneid23 <-convnafl_DOWN
geneid24 <-convf0_DOWN
geneid25 <-convf2_DOWN
geneid26 <-convf3_DOWN
geneid27 <-convf4_DOWN
geneid28 <-naflvf2_DOWN
geneid29 <-naflvf3_DOWN
geneid30 <-naflvf4_DOWN
geneid31 <-f0vf2_DOWN
geneid32 <-f0vf3_DOWN
geneid33 <-f0vf4_DOWN

geneid34 <-convnafl_DOWN_NEUT
geneid35 <-convf0_DOWN_NEUT
geneid36 <-convf2_DOWN_NEUT
geneid37 <-convf3_DOWN_NEUT
geneid38 <-convf4_DOWN_NEUT
geneid39 <-naflvf2_DOWN_NEUT
geneid40 <-naflvf3_DOWN_NEUT
geneid41 <-naflvf4_DOWN_NEUT
geneid42 <-f0vf2_DOWN_NEUT
geneid43 <-f0vf3_DOWN_NEUT
geneid44 <-f0vf4_DOWN_NEUT







                                        #Gene Ontology Analysis



                                        #Biological Processess


GO_BP <- enrichGO(gene = geneid11$gene,
                  keyType       = 'ENSEMBL',
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  minGSSize     = 5,
                  readable      = TRUE)

gobp1 <- as.data.frame(summary(GO_BP))
head(gobp1)

write.table(gobp1,"F0F1vF4_GO_BP_UP_0,01", sep="\t", quote=F)

plot1 <- dotplot(GO_BP, showCategory=15,title="BP: F0F1 vs F4 - Upregulated")
plot1

ggsave("F0F1vsF4_GO_BP_DOT_UP_0,01.pdf", plot1)


plot2 <- cnetplot(GO_BP, 5, colorEdge=T)
plot2

ggsave("F0F1vsF4_GO_BP_NET_UP_0,01.pdf", plot2)



                                        #Molecular Functions


GO_MF <- enrichGO(gene = geneid11$gene,
                keyType       = 'ENSEMBL',
                OrgDb         = org.Hs.eg.db,
                ont           = "MF",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                minGSSize     = 5,
                readable      = TRUE)

gomf1 <- as.data.frame(summary(GO_MF))
head(gobp1)

write.table(gomf1,"F0F1vsF4_GO_MF_UP_0,01", sep="\t", quote=F)

plot1 <- dotplot(GO_MF, showCategory=15,title="MF: F0F1 vs F4 - Upregulated")
plot1

ggsave("F0F1vsF4_GO_MF_DOT_UP_0,01.pdf", plot1)


plot2 <- cnetplot(GO_MF,showCategory= 5, colorEdge=T)
plot2

ggsave("F0F1vsNAFL_GO_MF_NET_UP_0,01.pdf", plot2)



                                        #Cellular Component



GO_CC <- enrichGO(gene = geneid11$gene,
                keyType       = 'ENSEMBL',
                OrgDb         = org.Hs.eg.db,
                ont           = "CC", ###Try to c("CC","BP","MF","KEGG")
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                minGSSize     = 5,
                readable      = TRUE)

gocc1 <- as.data.frame(summary(GO_CC))
head(gocc1)

colnames(gocc1)



write.table(gocc1,"F0F1vsF4_GO_CC_UP_0,01", sep="\t", quote=F)

plot1 <- dotplot(GO_CC, showCategory=5,title="CC: F0F1 vs F4 - Upregulated")
plot1

ggsave("F0F1vsF4_GO_CC_DOT_UP_0,01.pdf", plot1)


plot2 <- cnetplot(GO_CC, 5, colorEdge=T)
plot2

ggsave("ControlvsNAFL_GO_CC_NET_UP_0,01.pdf", plot2)



                                        #KEGG


kegg1 <- bitr(geneid11$gene, fromType = "ENSEMBL",
        toType = c("ENTREZID"),
        OrgDb = org.Hs.eg.db)



##names(kegg1) <- c("gene", "ENTREZID")
##keggfix <- merge(kegg1, test4, by="gene")
##keggfix <- keggfix[,-1]


kegg2 <- enrichKEGG(gene = kegg1$ENTREZID, ##Replace with keggfix in case its not working
                 organism     = 'hsa',
                 pvalueCutoff = 0.05,
                 qvalueCutoff  = 0.05,
                 minGSSize = 3)

head(kegg2)

write.table(kegg2,"F0F1vsF4_KEGG_UP_0,01", sep="\t", quote=F)

plot1 <- dotplot(kegg2, showCategory=15, title="KEGG: F0-F1 vs F4 - Upregulated")
plot1

ggsave("F0F1vsF4_KEGG_DOT_UP_0,01.pdf", plot1)

plot2 <- cnetplot(kegg2, 5, colorEdge=T)
plot2

ggsave("F0F1vsF4_KEGG_NET_UP_0,01.pdf", plot2)



## gProfiler GO/KEGG // Alternative method not use in Thesis

library(gprofiler2)

test <- gconvert(
geneid1$gene,
organism = "hsapiens",
target = "ENSG",
numeric_ns = "",
mthreshold = Inf,
filter_na = FALSE
)

head(test)

test1 <- gost(
geneid1$gene,   ###Try test after gconvert//geneid1 originally
organism = "hsapiens",
ordered_query = FALSE,
multi_query = FALSE,
significant = TRUE,
exclude_iea = FALSE,
measure_underrepresentation = FALSE,
evcodes = FALSE,
user_threshold = 0.05,
correction_method = c("g_SCS", "bonferroni", "fdr", "false_discovery_rate", "gSCS",
"analytical"),
domain_scope = c("annotated", "known", "custom", "custom_annotated"),
custom_bg = NULL,
numeric_ns = "",
sources = c("GO:BP","GO:MF","GO:CC"), #Add KEGG
as_short_link = FALSE)


plot1 <- dotplot(test1, showCategory=15,title="Control vs NAFL - Upregulated")
plot1

ggsave("CtrlvsNAFL_GO_KEGG_DOT_UP_0,01.pdf", plot1)


plot2 <- cnetplot(test1, 5, colorEdge=T)
plot2

ggsave("F0F1vsF4_GO_KEGG_NET_UP_0,01.pdf", plot2)








test1 <- gost(
xx2$gene,   ###Try test after gconvert//geneid1 originally
organism = "hsapiens",
ordered_query = FALSE,
multi_query = FALSE,
significant = TRUE,
exclude_iea = FALSE,
measure_underrepresentation = FALSE,
evcodes = FALSE,
user_threshold = 0.05,
correction_method = c("g_SCS", "bonferroni", "fdr", "false_discovery_rate", "gSCS",
"analytical"),
domain_scope = "annotated",
custom_bg = NULL,
numeric_ns = "",
sources = c("GO:BP","GO:MF","GO:CC"), #Add KEGG
as_short_link = FALSE)


head(test1$result, 5)


rownames(geneid1) <- NULL

print(geneid1, row.names=FALSE)

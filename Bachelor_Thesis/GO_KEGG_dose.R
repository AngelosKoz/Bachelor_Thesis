setwd('/home/aggelos/Desktop/Rotation/Files')


##library("PloGO2")

library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)

library(parallel)

sessionInfo()

##protcounts <- read.delim("protein", sep="\t", header=T)
##protcounts <- as.data.frame(t(protcounts))
##protcounts <- protcounts[-c(1:3),]

##patients <- read.delim("patorderstages", sep="\t", header=T)
##patNOf0f1 <- as.data.frame(patients[cc(62:95),-c(2:5)])
##names(patNOf0f1) <- "patient"
##patNOf0f1


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



geneid <- as.data.frame(convnafl[,1]) ##CHANGE when you input new set of genes
names(geneid) <- "gene"
geneid




test <- read.delim("MLproteindf", sep="\t", header=T) 
test1 <- test[-c(1:10,62:95), -c(18402:18405)] ##Remove patients not used 
test2 <- as.data.frame(t(test1))
test3 <- tibble::rownames_to_column(test2, "gene")
test4 <- merge(geneid, test3, by="gene")



                                        #Biological Processess


GO_BP <- enrichGO(gene = geneid$gene,              ###Change from geneid to test4 if different df
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

plot1 <- dotplot(GO_BP, showCategory=15)

plot2 <- cnetplot(GO_BP, 5)

dev.copy(plot1, "_GO_BP_DOT_0,01.pdf")
dev.off()

dev.copy(plot2, "_GO_BP_NET_0,01.pdf")
dev.off() 
write.table(gobp1,"_GO_BP_0,01", sep="\t", quote=F)

                                        #Molecular Functions


GO_MF <- enrichGO(gene = geneid$gene,
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

plot1 <- dotplot(GO_MF, showCategory=15)

plot2 <- cnetplot(GO_MF, 5)


dev.copy(plot1, "_GO_MF_DOT_0,01.pdf")
dev.off()

dev.copy(plot2, "_GO_MF_NET_0,01.pdf")
dev.off() 

write.table(gomf1,"_GO_MF_0,01", sep="\t", quote=F)

                                        #Cellular Component



GO_CC <- enrichGO(gene = geneid$gene,
                keyType       = 'ENSEMBL',
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                minGSSize     = 5,
                readable      = TRUE)

gocc1 <- as.data.frame(summary(GO_CC))
head(gocc1)

plot1 <- dotplot(GO_CC, showCategory=15)

plot2 <- cnetplot(GO_CC, 5)

dev.copy(plot1, "_GO_CC_DOT_0,01.pdf")

dev.off()

dev.copy(plot2, "_GO_CC_NET_0,01.pdf")
dev.off() 

write.table(gocc1,"_GO_CC_0,01", sep="\t", quote=F)

                                        #KEGG


kegg1 <- bitr(geneid$gene, fromType = "ENSEMBL",
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

plot1 <- dotplot(kegg2, showCategory=15)

plot2 <- cnetplot(kegg2, 5)


dev.copy(plot1, "_KEGG_DOT_0,01.pdf")
dev.off()

dev.copy(plot2, "_KEGG_NET_0,01.pdf")
dev.off()

write.table(kegg2,"_KEGG_0,01", sep="\t", quote=F)

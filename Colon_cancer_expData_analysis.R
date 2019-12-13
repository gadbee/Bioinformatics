rm(list=ls())
setwd('D:\\Projects\\colonCancer_expData_analysis')


# Read in TCGA_GTEx_TARGET expression data.###########################
library(data.table)
TcgaTargetGtex_expectedCount_log <- 
  fread('D:\\PublicData\\DataFromXena\\TcgaTargetGtex_gene_expected_count.gz',
        data.table=F, stringsAsFactors=F, sep='\t')
save(TcgaTargetGtex_expectedCount_log,
     file='D:\\Projects\\colonCancer_expData_analysis\\TcgaTargetGtex_expectedCount_log.RData')


# Extract colon sample ID from phenotype data.###################################################
library(stringr)
TcgaTargetGtex_phenotype <- 
  read.table('D:\\PublicData\\DatafromXena\\TcgaTargetGTEX_phenotype.txt.gz',
             header=T, stringsAsFactors=F, sep='\t')
colonTumor_TcgaGtex_phenotype <-
  TcgaTargetGtex_phenotype[which(TcgaTargetGtex_phenotype[,4]=='Colon'&
                                   str_detect(pattern='Primary Tumor',
                                              string=TcgaTargetGtex_phenotype[,5])),]
colonParaTumor_TcgaGtex_phenotype <-
  TcgaTargetGtex_phenotype[which(TcgaTargetGtex_phenotype[,4]=='Colon'&
                                   str_detect(pattern='Solid Tissue Normal',
                                              string=TcgaTargetGtex_phenotype[,5])),]
colonNormal_TcgaGtex_phenotype <- 
  TcgaTargetGtex_phenotype[which(TcgaTargetGtex_phenotype[,4]=='Colon'&
                                   str_detect(pattern='Normal Tissue',
                                              string=TcgaTargetGtex_phenotype[,5])),]
TcgaGtex_phenotype <- rbind(colonNormal_TcgaGtex_phenotype,
                                  colonParaTumor_TcgaGtex_phenotype,
                                  colonTumor_TcgaGtex_phenotype)
rm(colonTumor_TcgaGtex_phenotype, colonParaTumor_TcgaGtex_phenotype,
   colonNormal_TcgaGtex_phenotype, TcgaTargetGtex_phenotype)


# Extract colon expression data from TCGA_GTEx_TARGET expression data###############
TcgaGTEx_expectedCount_log <- 
  data.frame(EnsemblID=TcgaTargetGtex_expData_expectedCount_log$sample,
             TcgaTargetGtex_expData_expectedCount_log[, colnames(TcgaTargetGtex_expData_expectedCount_log)%in%
                                                        TcgaGtex_phenotype$sample])
#Attention! Attention!! R will transform the '-' in the colnames into '.'!!!
rm(TcgaTargetGtex_expData_expectedCount_log)


# Transform log(x+1) to expected count and filter low-expressed genes.#################
TcgaGTEx_expectedCount <- data.frame(EnsemblID=TcgaGTEx_expectedCount_log$EnsemblID,
                                     round(2^TcgaGTEx_expectedCount_log[,-1]-1),
                                     stringsAsFactors=F)
rm(TcgaGTEx_expectedCount_log)
TcgaGtex_expectedCount_filtered <- 
  TcgaGTEx_expectedCount[rowSums(TcgaGTEx_expectedCount[,-1])>1,]
rm(TcgaGTEx_expectedCount)


if (F) {
#DESeq analysis with EnsemblID
library(DESeq2)
library(BiocParallel)
dds <- DESeqDataSetFromMatrix(countData=TcgaGtex_expectedCount_filtered_matrix,
                              colData=groups, design=~condition)
register(SnowParam(6))
dds_out <- DESeq(dds, parallel=T)
save(dds_out, file='TcgaGtex_DEGs.RData')
TcgaGtex_DEGs_rslt <- results(dds_out)
write.csv(TcgaGtex_DEGs_rslt,
          file='TcgaGtex_DEGs_rslt.csv')
}


# Transform ensemblID to gene symbol, aggregate the dataframe by gene symbol and #########
#transform the data frame to matrix.
TcgaGtex_expectedCount_filtered$EnsemblID <- 
  unlist(strsplit(x=as.character(TcgaGtex_expectedCount_filtered$EnsemblID),
           split='\\.'))[seq(1,2*nrow(TcgaGtex_expectedCount_filtered),2)]
ensemblID_gene <- read.table('D:\\PublicData\\DataFromXena\\gencode.v23.annotation.gene.probemap',
                             header=T, stringsAsFactors=F)
ensemblID_gene <- ensemblID_gene[,c(1,2)]
ensemblID_gene$id <- unlist(strsplit(x=ensemblID_gene$id, 
                                     split='\\.'))[seq(1, 2*nrow(ensemblID_gene),2)]
names(ensemblID_gene)[1] <- 'EnsemblID'
TcgaGtex_expectedCount_filtered_gene <-
  merge(TcgaGtex_expectedCount_filtered, ensemblID_gene, by='EnsemblID')
rm(TcgaGtex_expectedCount_filtered, ensemblID_gene)
TcgaGtex_expectedCount_filtered_gene <-
  TcgaGtex_expectedCount_filtered_gene[,-1]
TcgaGtex_expectedCount_filtered_gene_agr <-
  aggregate(TcgaGtex_expectedCount_filtered_gene[,-length(TcgaGtex_expectedCount_filtered_gene)],
            by=list(TcgaGtex_expectedCount_filtered_gene$gene), FUN=sum)
save(TcgaGtex_expectedCount_filtered_gene_agr,
     file='TcgaGtex_expectedCount_filtered_gene_agr.RData')
rm(TcgaGtex_expectedCount_filtered_gene)
rownames(TcgaGtex_expectedCount_filtered_gene_agr) <-
  TcgaGtex_expectedCount_filtered_gene_agr$Group.1
TcgaGtex_expectedCount_filtered_gene_agr <-
  TcgaGtex_expectedCount_filtered_gene_agr[,-1]
TcgaGtex_expectedCount_filtered_gene_agr_matrix <- 
  as.matrix(TcgaGtex_expectedCount_filtered_gene_agr)


# 3 groups analysis #######################################################################

library(DESeq2)
library(BiocParallel)
register(SnowParam(12))

#Establish grouping vector and grouping matrix with 3 groups.
TcgaGtex_phenotype$sample <- gsub(pattern='-', replacement='.',
                                        x=TcgaGtex_phenotype$sample)
TcgaGtex_phenotype <-
  TcgaGtex_phenotype[TcgaGtex_phenotype$sample%in%
                       colnames(TcgaGtex_expectedCount_filtered_gene_agr),
                     c(1,5)]
names(TcgaGtex_phenotype)[2] <- 'condition'
TcgaGtex_phenotype_with_3groups <- merge(data.frame(sample=colnames(TcgaGtex_expectedCount_filtered_gene_agr)),
                                         TcgaGtex_phenotype, by='sample', sort=F)
rm(TcgaGtex_phenotype)

#PCA analysis with 3 groups.
gene_dds_3groups <- DESeqDataSetFromMatrix(countData=TcgaGtex_expectedCount_filtered_gene_agr_matrix,
                                           colData=TcgaGtex_phenotype_with_3groups,
                                           design=~condition)
gene_dds_out_3groups <- DESeq(gene_dds_3groups, parallel=T)
save(gene_dds_out_3groups, file='TcgaGtex_gene_DEGs_3groups.RData')

vsd_3groups <- vst(gene_dds_out_3groups, blind=F)
plotPCA(vsd_3groups, 'condition')

#Solid Tissue Normal vs Normal Tissue
TcgaGtex_phenotyp_ParatumorVsNormal <- 
  subset(TcgaGtex_phenotype_with_3groups,
         condition=='Solid Tissue Normal'|condition=='Normal Tissue')
TcgaGtex_expectedCount_filtered_gene_agr_ParatumorVsNormal <-
  TcgaGtex_expectedCount_filtered_gene_agr[,colnames(TcgaGtex_expectedCount_filtered_gene_agr)
                                                     %in%TcgaGtex_phenotyp_ParatumorVsNormal$sample]
TcgaGtex_expectedCount_filtered_gene_agr_ParatumorVsNormal_matrix <-
  as.matrix(TcgaGtex_expectedCount_filtered_gene_agr_ParatumorVsNormal)
gene_dds_ParatumorVsNormal <- DESeqDataSetFromMatrix(countData=TcgaGtex_expectedCount_filtered_gene_agr_ParatumorVsNormal_matrix,
                                                     colData=TcgaGtex_phenotyp_ParatumorVsNormal,
                                                     design=~condition)
gene_dds_out_ParatumorVsNormal <- DESeq(gene_dds_ParatumorVsNormal, parallel=T)
save(gene_dds_out_ParatumorVsNormal, file='gene_dds_out_ParatumorVsNormal.RData')
result_ParatumorVsNormal <- results(gene_dds_out_ParatumorVsNormal,
                                    contrast=c('condition', 'Solid Tissue Normal',
                                               'Normal Tissue'))
write.csv(result_ParatumorVsNormal, file='result_ParatumorVsNormal.csv')
vsd_ParatumorVsNormal <- vst(gene_dds_out_ParatumorVsNormal, blind=F)
plotPCA(vsd_ParatumorVsNormal,'condition')

#Primary Tumor vs Normal Tissue
TcgaGtex_phenotyp_TumorVsNormal <- subset(TcgaGtex_phenotype_with_3groups,
                                              condition=='Primary Tumor'|
                                                condition=='Normal Tissue')
TcgaGtex_expectedCount_filtered_gene_agr_TumorVsNormal <-
  TcgaGtex_expectedCount_filtered_gene_agr[,colnames(TcgaGtex_expectedCount_filtered_gene_agr)
                                           %in%TcgaGtex_phenotyp_TumorVsNormal$sample]
TcgaGtex_expectedCount_filtered_gene_agr_TumorVsNormal_matrix <-
  as.matrix(TcgaGtex_expectedCount_filtered_gene_agr_TumorVsNormal)
gene_dds_TumorVsNormal <- DESeqDataSetFromMatrix(countData=TcgaGtex_expectedCount_filtered_gene_agr_TumorVsNormal_matrix,
                                                     colData=TcgaGtex_phenotyp_TumorVsNormal,
                                                     design=~condition)
gene_dds_out_TumorVsNormal <- DESeq(gene_dds_TumorVsNormal, parallel=T)
save(gene_dds_out_TumorVsNormal, file='gene_dds_out_TumorVsNormal.RData')
result_TumorVsNormal <- results(gene_dds_out_TumorVsNormal,
                                contrast=c('condition', 'Primary Tumor', 'Normal Tissue'))
write.csv(result_TumorVsNormal, file='result_TumorVsNormal.csv')
vsd_TumorVsNormal <- vst(gene_dds_out_TumorVsNormal, blind=F)
plotPCA(vsd_TumorVsNormal, 'condition')

#Primary Tumor vs Solid Tissue Normal
TcgaGtex_phenotype_TumorVsParatumor <-
  TcgaGtex_phenotype_with_3groups[TcgaGtex_phenotype_with_3groups$
                                          condition=='Primary Tumor'
                                          |TcgaGtex_phenotype_with_3groups$
                                          condition=='Solid Tissue Normal',]
TcgaGtex_expectedCount_filtered_gene_agr_TumorVsParatumor <-
  TcgaGtex_expectedCount_filtered_gene_agr[,colnames(TcgaGtex_expectedCount_filtered_gene_agr)
                                                         %in%TcgaGtex_phenotype_TumorVsParatumor$sample]
TcgaGtex_expectedCount_filtered_gene_agr_TumorVsParatumor_matrix <-
  as.matrix(TcgaGtex_expectedCount_filtered_gene_agr_TumorVsParatumor)
gene_dds_TumorVsParatumor <- 
  DESeqDataSetFromMatrix(TcgaGtex_expectedCount_filtered_gene_agr_TumorVsParatumor_matrix,
                         colData=TcgaGtex_phenotype_TumorVsParatumor,
                         design=~condition)
gene_dds_out_TumorVsParatumor <- DESeq(gene_dds_TumorVsParatumor, parallel=T)
save(gene_dds_out_TumorVsParatumor, file='gene_dds_out_TumorVsParatumor.RData')
result_TumorVsParatumor <- results(gene_dds_out_TumorVsParatumor,
                                   contrast=c('condition','Primary Tumor', 'Solid Tissue Normal'))
write.csv(result_TumorVsParatumor, file='result_TumorVsParatumor.csv')
vsd_TumorVsParatumor <- vst(gene_dds_out_TumorVsParatumor, blind=F)
plotPCA(vsd_TumorVsParatumor, 'condition')

# 3 groups analysis end 


# Enrichment analysis ###################################################################
library(clusterProfiler)
library(org.Hs.eg.db)

getGeneList4GoKeggEnrichment <- function(resultDataFile){
  #Transform the DEG data to the data for GO and KEGG enrichment analysis.
  result <- read.csv(file=resultDataFile, header=T, stringsAsFactors=F)
  names(result)[1] <- 'GeneSymbol'
  DEGs <- subset(x=result, abs(log2FoldChange)>2&padj<0.05)$GeneSymbol
  symbol_ensemblID <- read.table('D:\\PublicData\\DataFromXena\\gencode.v23.annotation.gene.probemap',
                                 header=T, stringsAsFactors=F)
  DE_ensemblID <- merge(data.frame(gene=DEGs),symbol_ensemblID, by.x='gene')
  DE_ensemblID <- unique(DE_ensemblID$id)
  DE_ensemblID <- unlist(strsplit(DE_ensemblID, split='\\.'))[seq(1,2*length(DE_ensemblID),2)]
  return(DE_ensemblID)
}

getGeneList4GseaEnrichment <- function(filePath){
  #Transform the DEG data to the data for GSEA enrichment analysis.
  resultFile <- read.table(filePath, header=T, stringsAsFactors=F, sep=',')
  forGseaData <- resultFile[,c(1,3)]
  names(forGseaData)[1] <- 'SYMBOL'
  EntrezID <- bitr(geneID=forGseaData[,1], fromType='SYMBOL', 
                   toType='ENTREZID', OrgDb=org.Hs.eg.db)
  logFC_info <- merge(forGseaData, EntrezID, by='SYMBOL')
  logFC <- logFC_info[,2]
  names(logFC) <- logFC_info[,3]
  logFC_sorted <- sort(logFC, decreasing=T)
  return(logFC_sorted)
}

GO_KEGG_enrichment <- function(geneset, dataname, pvalueCutoff, qvalueCutoff){
  ego_cc <- enrichGO(gene=geneset,
                     keyType='ENSEMBL',
                     OrgDb=org.Hs.eg.db,
                     ont='CC',
                     pAdjustMethod='BH',
                     pvalueCutoff=pvalueCutoff,
                     qvalueCutoff=qvalueCutoff,
                     readable=T)
  
  ego_bp <- enrichGO(gene=geneset,
                     keyType='ENSEMBL',
                     OrgDb=org.Hs.eg.db,
                     ont='BP',
                     pAdjustMethod='BH',
                     pvalueCutoff=pvalueCutoff,
                     qvalueCutoff=qvalueCutoff,
                     readable=T)
  
  ego_mf <- enrichGO(gene=geneset,
                     keyType='ENSEMBL',
                     OrgDb=org.Hs.eg.db,
                     ont='MF',
                     pAdjustMethod='BH',
                     pvalueCutoff=pvalueCutoff,
                     qvalueCutoff=qvalueCutoff,
                     readable=T)
  
  entrezID <- bitr(geneset,
                   fromType='ENSEMBL',
                   toType='ENTREZID',
                   OrgDb=org.Hs.eg.db)
  
  kk <- enrichKEGG(gene=entrezID$ENTREZID,
                   organism='human',
                   pvalueCutoff=pvalueCutoff,
                   qvalueCutoff=qvalueCutoff)
  
 
  results <- list()
  results$CC <- ego_cc
  results$BP <- ego_bp
  results$MF <- ego_mf
  results$KK <- kk
  
  return(results)
}

GseaEnrichment <- function(geneSet, dataname, nPerm, pvalueCutoff){
  gsea_cc <- gseGO(geneList=geneSet,
                   OrgDb=org.Hs.eg.db,
                   ont="CC",
                   nPerm=nPerm,
                   maxGSSize=500,
                   pvalueCutoff=pvalueCutoff,
                   verbose=FALSE)
  
  gsea_bp <- gseGO(geneList=geneSet,
                   OrgDb=org.Hs.eg.db,
                   ont="BP",
                   nPerm=nPerm,
                   maxGSSize=500,
                   pvalueCutoff=pvalueCutoff,
                   verbose=FALSE)
  
  gsea_mf <- gseGO(geneList=geneSet,
                   OrgDb=org.Hs.eg.db,
                   ont="MF",
                   nPerm=nPerm,
                   maxGSSize=500,
                   pvalueCutoff=pvalueCutoff,
                   verbose=FALSE)
  
  gsea_kegg <- gseKEGG(geneList=geneSet,
                       organism='hsa',
                       nPerm=nPerm,
                       minGSSize=15,
                       pvalueCutoff=pvalueCutoff,
                       verbose=FALSE)

  results <- list()
  results$CC <- gsea_cc
  results$BP <- gsea_bp
  results$MF <- gsea_mf
  results$KK <- gsea_kegg
  return(results)
}

# Enrichment of TumorVsNormal DE genes
TumorVsNormal_DE_genes <- getGeneList4GoKeggEnrichment('result_TumorVsNormal.csv')
GO_KEGG_TumorVsNormal <-  GO_KEGG_enrichment(TumorVsNormal_DE_genes,'TumorVsNormal', 0.01, 0.01)
gsea_TumorVsNormal_genes <- getGeneList4GseaEnrichment('result_TumorVsNormal.csv')
GSEA_TumorVsNormal <- GseaEnrichment(gsea_TumorVsNormal_genes, 'TumorVsNormal', 1000, 0.05)

# Enrichment of TumorVsParatumor DE genes
TumorVsParatumor_DE_genes <- getGeneList4GoKeggEnrichment('result_TumorVsParatumor.csv')
GO_KEGG_TumorVsParatumor <- GO_KEGG_enrichment(TumorVsParatumor_DE_genes,
                                               'TumorVsParatumor', 0.01, 0.01)
gsea_TumorVsParatumor_genes <- getGeneList4GseaEnrichment('result_TumorVsParatumor.csv')
gsea_TumorVsParatumor <- GseaEnrichment(gsea_TumorVsParatumor_genes, 'TumorVsParatumor', 1000, 0.05)

# Enrichment of ParatumorVsNormal DE genes
ParatumorVsNormal_DE_genes <- getGeneList4GoKeggEnrichment('result_ParatumorVsNormal.csv')
GO_KEGG_ParatumorVsNormal <- GO_KEGG_enrichment(ParatumorVsNormal_DE_genes,
                                                'ParatumorVsNormal', 0.01, 0.01)
gsea_ParatumorVsNormal_genes <- getGeneList4GseaEnrichment('result_ParatumorVsNormal.csv')
gsea_ParatumorVsNormal <- GseaEnrichment(gsea_ParatumorVsNormal_genes, 'ParatumorVsNormal', 1000, 0.05)

# Homodromous changed genes enrichmeng
result_ParatumorVsNormal <- read.csv('result_ParatumorVsNormal.csv', header=T, stringsAsFactors=F)
result_TumorVsParatumor <- read.csv('result_TumorVsParatumor.csv', header=T, stringsAsFactors=F)
merged_result <- merge(result_ParatumorVsNormal, result_TumorVsParatumor, by='X')
merged_logFC_result <- data.frame(gene=merged_result$X,
                                  logFC.PN=merged_result$log2FoldChange.x,
                                  logFC.TP=merged_result$log2FoldChange.y)
homodromous_changed_genes_info <- subset(merged_logFC_result, (logFC.PN>0&logFC.TP>0)|
                                      (logFC.PN<0&logFC.TP<0))
result_TumorVsNormal <- read.csv('result_TumorVsNormal.csv', header=T, stringsAsFactors=F)
homodromous_changed_genes_TN <- subset(result_TumorVsNormal,
                                       result_TumorVsNormal$X %in% homodromous_changed_genes_info$gene)
homodromous_changed_genes_TN <- data.frame(gene=homodromous_changed_genes_TN$X, stringsAsFactors=F,
                                           logFC=homodromous_changed_genes_TN$log2FoldChange)
homodromous_changed_genes_ENTREZID <- bitr(homodromous_changed_genes_TN$gene, fromType='SYMBOL',
                                           toType='ENTREZID', OrgDb=org.Hs.eg.db)
names(homodromous_changed_genes_TN)[1] <- 'SYMBOL'
homodromous_changed_genes_TN <- merge(homodromous_changed_genes_ENTREZID, homodromous_changed_genes_TN,
                                        by='SYMBOL')
logFC_TN <- homodromous_changed_genes_TN$logFC
names(logFC_TN) <- homodromous_changed_genes_TN$ENTREZID
logFC_TN <- sort(logFC_TN, decreasing=T)

gsea_homodrmous_changed_genes <- GseaEnrichment(logFC_TN, 'homodrmous_changed_genes', 1000, 0.05)

ParatumorVsNormal_DE_genes <- subset(result_ParatumorVsNormal, padj<0.05&abs(log2FoldChange)>2)
TumorVsParatumor_DE_genes <- subset(result_TumorVsParatumor, padj<0.05&abs(log2FoldChange)>2)
temp <- merge(ParatumorVsNormal_DE_genes, TumorVsParatumor_DE_genes, by='X', sort=T)
homodromous_changed_sig_genes <- subset(temp, (log2FoldChange.x<0&log2FoldChange.y<0)|
                                      (log2FoldChange.x>0&log2FoldChange.y>0))
symbol_ensemblID <- read.table('D:\\PublicData\\DataFromXena\\gencode.v23.annotation.gene.probemap',
                               header=T, stringsAsFactors=F)
HCG_ensemblID_info <- merge(data.frame(gene=homodromous_changed_sig_genes[,1], stringsAsFactors=F),
                       symbol_ensemblID, by.x='gene')
HCG_ensemblID <- unique(HCG_ensemblID_info$id)
HCG_ensemblID <- unlist(strsplit(HCG_ensemblID, split='\\.'))[seq(1,2*length(HCG_ensemblID),2)]
GO_KEGG_HCG <- GO_KEGG_enrichment(HCG_ensemblID, 'HCG', 0.01, 0.01)

# Enrichment analysis based on gmt file
library(magrittr)
wpgmtfile <- system.file("extdata/wikipathways-20180810-gmt-Homo_sapiens.gmt", package="clusterProfiler")
wp2gene <- read.gmt(wpgmtfile)
wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME

ewp_TumorVsNormal <- enricher(names(gsea_TumorVsNormal_genes), 
                              TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
ewp_TumorVsParatumor <- enricher(names(gsea_TumorVsParatumor_genes),
                                 TERM2GENE=wpid2gene, TERM2NAME=wpid2name)
ewp_ParatumorVsNormal <- enricher(names(gsea_ParatumorVsNormal_genes),
                                 TERM2GENE=wpid2gene, TERM2NAME=wpid2name)

gsea_gmt_TumorVsNormal <- GSEA(gsea_TumorVsNormal_genes,
                               TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
gsea_gmt_TumorVsParatumor <- GSEA(gsea_TumorVsParatumor_genes,
                                  TERM2GENE=wpid2gene, TERM2NAME=wpid2name)
gsea_gmt_ParatumorVsNormal <- GSEA(gsea_ParatumorVsNormal_genes,
                                   TERM2GENE=wpid2gene, TERM2NAME=wpid2name)

# Enrichment analysis end 

# Draw enrichment plot ###########################################################
batch_dotplot <- function(group){
  group_list <- get(group)
  pdf(file=paste(group, '.pdf', sep=''))
  for (i in 1:4) {
    p <- dotplot(group_list[[i]], title=paste(group, names(group_list)[i]))
    print(p)
  }
  dev.off()
}

batch_dotplot('GO_KEGG_TumorVsNormal')
batch_dotplot('GSEA_TumorVsNormal')
batch_dotplot('GO_KEGG_TumorVsParatumor')
batch_dotplot('gsea_TumorVsParatumor')
batch_dotplot('GO_KEGG_ParatumorVsNormal')
batch_dotplot('gsea_ParatumorVsNormal')
batch_dotplot('GO_KEGG_HCG')
batch_dotplot('gsea_homodrmous_changed_genes')

pdf('enricher.pdf')
dotplot(ewp_TumorVsParatumor, title='ewp_TumorVsParatumor', showCategory=20)
dotplot(ewp_ParatumorVsNormal, title='ewp_ParatumorVsNormal', showCategory=20)
dotplot(ewp_TumorVsNormal, title='ewp_TumorVsNormal', showCategory=20)
dev.off()

library(enrichplot)
pdf(file='gsea_enrichment.pdf')
gseaplot2(gsea_gmt_TumorVsParatumor, geneSetID=1, title='TumorVsParatumor')
gseaplot2(gsea_gmt_ParatumorVsNormal, geneSetID=1, title='ParatumorVsNormal')
gseaplot2(gsea_gmt_TumorVsNormal, geneSetID=1, title='TumorVsNormal')
dev.off()

#Draw volcano plot
resultProcessing <- function(result, sigLevel, logFC){
  log2_padj <- data.frame(geneSymbol=rownames(result),
                          log2FC=result$log2FoldChange,
                          padj=result$padj)
  log2_padj <- na.omit(log2_padj)
  log2_padj$sig <- 1
  log2_padj$sig[log2_padj$padj<sigLevel&log2_padj$log2FC>=logFC] <- 2
  log2_padj$sig[log2_padj$padj<sigLevel&log2_padj$log2FC<=-logFC] <- 3
  log2_padj$showName <- F
  log2_padj$showName[which(log2_padj$padj<sigLevel)] <- T
  return(log2_padj)
}

volcanoPlot <-  function(pdfName,log2_padj, sigThreshold, logFCthreshold){
  pdf(pdfName)
  plot(log2_padj$log2FC, -log10(log2_padj$padj),
       xlab='log2FC', ylab='-log10(padj)',
       col=log2_padj$sig, xlim=c(-10,10), ylim=c(0,200))
  text(log2_padj[which(log2_padj$sig!=1),2],
       -log10(log2_padj[which(log2_padj$sig!=1),3]),
       log2_padj[which(log2_padj$sig!=1),1],
       col=log2_padj[which(log2_padj$sig!=1),4],
       pos=3, cex=0.8)
  abline(v=c(-logFCthreshold,logFCthreshold), h=-log10(sigThreshold), col='grey', lty=2)
  dev.off()
}
a <- resultProcessing(result_ParatumorVsNormal, 0.001, 5)
b <- resultProcessing(result_TumorVsNormal, 0.05, 2)
c <- resultProcessing(result_TumorVsParatumor, 0.05, 2)
volcanoPlot('ParatumorVsNormal.pdf', a, 0.001, 5)
volcanoPlot('TumorVsNormal.pdf', b, 0.05, 2)
volcanoPlot('TumorVsParatumor.pdf', c, 0.05, 2)
# Draw enrichment plot end


# Read in GDC colon cancer expression data ##############################################
library(data.table)
GDC_colonCancer_ExpData_log <- fread('D:\\PublicData\\DataFromXena\\TCGA-COAD.htseq_counts.tsv',
                                     data.table=F, stringsAsFactors=F, sep='\t')
GDC_colonCancer_ExpData <- data.frame(Ensembl_ID=GDC_colonCancer_ExpData_log$Ensembl_ID,
                                      round(2^GDC_colonCancer_ExpData_log[,-1]-1),
                                      stringsAsFactors=F)
rm(GDC_colonCancer_ExpData_log)
library(stringr)
GDC_colonCancer_ExpData$Ensembl_ID <- unlist(strsplit(x=GDC_colonCancer_ExpData$Ensembl_ID,
                                                      split='\\.'))[seq(1, length(GDC_colonCancer_ExpData$Ensembl_ID)*2,2)]

ensemblID_gene_GDC <- read.table('D:\\PublicData\\DataFromXena\\gencode.v23.annotation.gene.probemap',
                                 header=T, stringsAsFactors=F)
ensemblID_gene_GDC$id <- unlist(strsplit(x=ensemblID_gene_GDC$id,
                                         split='\\.'))[seq(1, length(ensemblID_gene_GDC$id)*2, 2)]
ensemblID_gene_GDC <- data.frame(Ensembl_ID=ensemblID_gene_GDC$id,
                                 GeneSymbol=ensemblID_gene_GDC$gene,
                                 stringsAsFactors=F)
GDC_colonCancer_ExpData_gene <- merge(GDC_colonCancer_ExpData, ensemblID_gene_GDC,
                                      by='Ensembl_ID', all.x=T)
GDC_colonCancer_ExpData_gene <- data.frame(GeneSymbol=GDC_colonCancer_ExpData_gene$GeneSymbol,
                                           GDC_colonCancer_ExpData_gene[,c(-1,-length(GDC_colonCancer_ExpData_gene))],
                                           stringsAsFactors=F)
GDC_colonCancer_ExpData_gene_agr <- aggregate(GDC_colonCancer_ExpData_gene[,-1],
                                              by=list(GDC_colonCancer_ExpData_gene$GeneSymbol),
                                              FUN=sum)
GDC_colonCancer_ExpData_gene_agr_matrix <- as.matrix(x=GDC_colonCancer_ExpData_gene_agr)
row.names(GDC_colonCancer_ExpData_gene_agr_matrix) <- GDC_colonCancer_ExpData_gene_agr_matrix[,1]
GDC_colonCancer_ExpData_gene_agr_matrix <- GDC_colonCancer_ExpData_gene_agr_matrix[,-1]
GDC_colonCancer_ExpData_gene_agr_matrix <- apply(GDC_colonCancer_ExpData_gene_agr_matrix,
                                                 2, as.numeric)
# Above step will eliminate rownames of matrix. SHIT!
row.names(GDC_colonCancer_ExpData_gene_agr_matrix) <- GDC_colonCancer_ExpData_gene_agr$Group.1

cli_info <- read.csv('D:\\PublicData\\DataFromXena\\GDC_colon_clinical_info.csv',header=T,
                       stringsAsFactors=F)
stage <- data.frame(sample_id=cli_info[,1], stage=cli_info$tumor_stage.diagnoses,
                    stringsAsFactors=F)
stage$sample_id <- gsub(x=stage$sample_id, pattern='-', replacement='\\.')
sample_id <- data.frame(sample_id=names(GDC_colonCancer_ExpData_gene_agr)[-1], stringsAsFactors=F)
group <- merge(sample_id, stage, by='sample_id', all.x=T, sort=F)
group <- group[group$stage!=''&group$stage!='not reported',]
GDC_colonCancer_ExpData_gene_agr_matrix <-
  GDC_colonCancer_ExpData_gene_agr_matrix[,which(colnames(GDC_colonCancer_ExpData_gene_agr_matrix)%in%group$sample_id)]

save(GDC_colonCancer_ExpData_gene_agr_matrix, file='GDC_colonCancer_ExpData_gene_agr_matrix.RData')
group$stage <- gsub(x=group$stage, pattern='stage ia', replacement='stage i')
group$stage <- gsub(x=group$stage, pattern='stage iia', replacement='stage ii')
group$stage <- gsub(x=group$stage, pattern='stage iib', replacement='stage ii')
group$stage <- gsub(x=group$stage, pattern='stage iic', replacement='stage ii')
group$stage <- gsub(x=group$stage, pattern='stage iiia', replacement='stage iii')
group$stage <- gsub(x=group$stage, pattern='stage iiib', replacement='stage iii')
group$stage <- gsub(x=group$stage, pattern='stage iiic', replacement='stage iii')
group$stage <- gsub(x=group$stage, pattern='stage iva', replacement='stage iv')
group$stage <- gsub(x=group$stage, pattern='stage ivb', replacement='stage iv')

library(DESeq2)
library(BiocParallel)
register(SnowParam(12))
dds <- DESeqDataSetFromMatrix(countData=GDC_colonCancer_ExpData_gene_agr_matrix,
                              colData=group, design=~stage)
dds_out <- DESeq(dds, parallel=T,)
save(dds_out, file='GDC_colonCancer_stage_dds_out.RData')
vsd <- vst(dds_out, blind=F)
plotPCA(vsd, 'stage')

Stage2VsStage1 <- results(object=dds_out, contrast=c('stage', 'stage ii', 'stage i'))
Stage3VsStage2 <- results(object=dds_out, contrast=c('stage', 'stage iii', 'stage ii'))
stage4VsStage3 <- results(object=dds_out, contrast=c('stage', 'stage iv', 'stage iii'))
tmp <- merge(Stage2VsStage1, Stage3VsStage2, by='row.names')
row.names(tmp) <- tmp$Row.names
tmp <- tmp[,-1]
stage_logFC <- merge(tmp, stage4VsStage3, by='row.names')
stage_logFC <- data.frame(GeneSymbol=stage_logFC$Row.names,
                          logFC2Vs1=stage_logFC$log2FoldChange.x,
                          logFC3Vs2=stage_logFC$log2FoldChange.y,
                          logFC4Vs3=stage_logFC$log2FoldChange)
stage_logFC_cochange <- stage_logFC[(stage_logFC$logFC2Vs1>0.5&
                                       stage_logFC$logFC3Vs2>0.5&
                                       stage_logFC$logFC4Vs3>0.5)|
                                      (stage_logFC$logFC2Vs1<(-0.5)&
                                         stage_logFC$logFC3Vs2<(-0.5)&
                                         stage_logFC$logFC4Vs3<(-0.5)),]
stage_logFC_cochange <- na.omit(stage_logFC_cochange)
stage_logFC_cochange$mean <- apply(stage_logFC_cochange[,2:4],1,mean)


cochange_ensemblID <- merge(stage_logFC_cochange, ensemblID_gene_GDC, by='GeneSymbol')
cochange_ensemblID <- cochange_ensemblID$Ensembl_ID
GO_KEGG_cochange <-  GO_KEGG_enrichment(cochange_ensemblID,'cochange', 0.05, 0.05)

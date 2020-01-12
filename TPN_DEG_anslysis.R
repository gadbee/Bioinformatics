rm(list=ls())
setwd('D:\\Projects\\colonCancer_expData_analysis')
load('D:\\Projects\\colonCancer_expData_analysis\\TcgaTargetGtex_expectedCount_log.RData')

# Extract colon sample ID from phenotype data. ###################################################
library(stringr)
TcgaTargetGtex_phenotype <- 
  read.table('D:\\PublicData\\DatafromXena\\TcgaTargetGTEX_phenotype.txt.gz',
             header=T, stringsAsFactors=F, sep='\t')
colonTumor_TcgaGtex_phenotype <-
  TcgaTargetGtex_phenotype[which(TcgaTargetGtex_phenotype[,4]=='Colon'& 
                                   str_detect(pattern='Primary Tumor', string=TcgaTargetGtex_phenotype[,5])),]
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
#Attention! Attention!! R will transform the '-' in the colnames into '.' automatically!!!
rm(TcgaTargetGtex_expData_expectedCount_log)

# Transform log(x+1) to expected count and filter low-expressed genes.#################
TcgaGTEx_expectedCount <- data.frame(EnsemblID=TcgaGTEx_expectedCount_log$EnsemblID,
                                     round(2^TcgaGTEx_expectedCount_log[,-1]-1),
                                     stringsAsFactors=F)
rm(TcgaGTEx_expectedCount_log)
TcgaGtex_expectedCount_filtered <- 
  TcgaGTEx_expectedCount[rowSums(TcgaGTEx_expectedCount[,-1])>1,]
rm(TcgaGTEx_expectedCount)

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
names(ensemblID_gene) <- c('EnsemblID', 'GeneSymbol')
TcgaGtex_expectedCount_filtered_gene <-
  merge(TcgaGtex_expectedCount_filtered, ensemblID_gene, by='EnsemblID')
rm(TcgaGtex_expectedCount_filtered)
TcgaGtex_expectedCount_filtered_gene <-
  TcgaGtex_expectedCount_filtered_gene[,-1]
TcgaGtex_expectedCount_filtered_gene_agr <-
  aggregate(TcgaGtex_expectedCount_filtered_gene[,-length(TcgaGtex_expectedCount_filtered_gene)],
            by=list(TcgaGtex_expectedCount_filtered_gene$GeneSymbol), FUN=sum)
save(TcgaGtex_expectedCount_filtered_gene_agr,
     file='TcgaGtex_expectedCount_filtered_gene_agr.RData')
rm(TcgaGtex_expectedCount_filtered_gene)
rownames(TcgaGtex_expectedCount_filtered_gene_agr) <- TcgaGtex_expectedCount_filtered_gene_agr$Group.1
TcgaGtex_expectedCount_filtered_gene_agr <- TcgaGtex_expectedCount_filtered_gene_agr[,-1]
TcgaGtex_expectedCount_filtered_gene_agr_matrix <- as.matrix(TcgaGtex_expectedCount_filtered_gene_agr[,-1])

# DEG analysis #######################################################################

library(DESeq2)
library(BiocParallel)
register(SnowParam(12))

#Establish grouping vector and grouping matrix with 3 groups.
TcgaGtex_phenotype$sample <- gsub(pattern='-', replacement='.',
                                  x=TcgaGtex_phenotype$sample)
TcgaGtex_phenotype <- TcgaGtex_phenotype[TcgaGtex_phenotype$sample%in%
                                           colnames(TcgaGtex_expectedCount_filtered_gene_agr), c(1,5)]
names(TcgaGtex_phenotype)[2] <- 'condition'
TcgaGtex_phenotype_with_3groups <- TcgaGtex_phenotype[match(colnames(TcgaGtex_expectedCount_filtered_gene_agr_matrix), TcgaGtex_phenotype$sample),]
rm(TcgaGtex_phenotype)

#DESeq analysis and PCA.
gene_dds_3groups <- DESeqDataSetFromMatrix(countData=TcgaGtex_expectedCount_filtered_gene_agr_matrix,
                                           colData=TcgaGtex_phenotype_with_3groups,
                                           design=~condition)
gene_dds_out_3groups <- DESeq(gene_dds_3groups, parallel=T)
save(gene_dds_out_3groups, file='TcgaGtex_gene_DEGs_3groups.RData')
result_ParatumorVsNormal <- results(gene_dds_out_3groups,
                                    contrast=c('condition', 'Solid Tissue Normal',
                                               'Normal Tissue'))
result_TumorVsParatumor <- results(gene_dds_out_3groups,
                                   contrast=c('condition','Primary Tumor', 'Solid Tissue Normal'))
vsd_3groups <- vst(gene_dds_out_3groups, blind=F)
a <- plotPCA(vsd_3groups, 'condition')

# Remove outliers.
PCAdata <- a$data
outlier1 <- PCAdata[PCAdata$PC1>75&PCAdata$PC2<(-60), 5]
outlier2 <- PCAdata[PCAdata$PC1>0&PCAdata$group=='Solid Tissue Normal', 5]
outliers <- c(as.character(outlier1), as.character(outlier2))
TcgaGtex_expectedCount_filtered_gene_agr_outlierRemoved <- 
  TcgaGtex_expectedCount_filtered_gene_agr[,!colnames(TcgaGtex_expectedCount_filtered_gene_agr)%in%outliers]
TcgaGtex_expectedCount_filtered_gene_agr_outlierRemoved_matrix <- 
  as.matrix(TcgaGtex_expectedCount_filtered_gene_agr_outlierRemoved[,-1])
row.names(TcgaGtex_expectedCount_filtered_gene_agr_outlierRemoved_matrix) <-
  row.names(TcgaGtex_expectedCount_filtered_gene_agr_outlierRemoved)
TcgaGtex_phenotype_with_3groups_outlierRemoved <-
  TcgaGtex_phenotype_with_3groups[!TcgaGtex_phenotype_with_3groups$sample%in%outliers,]

# DESeq analysis of outlier removed data.
library(DESeq2)
library(parallel)
dds <- DESeqDataSetFromMatrix(countData=TcgaGtex_expectedCount_filtered_gene_agr_outlierRemoved_matrix,
                              colData=TcgaGtex_phenotype_with_3groups_outlierRemoved,
                              design=~condition)
dds_out <- DESeq(dds, parallel=T)
save(dds_out, file='TcgaGtex_gene_DEGs_3groups_outlierRemoved.RData')
result_ParatumorVsNormal_OR <- results(dds_out, contrast=c('condition',
                                                           'Solid Tissue Normal',
                                                           'Normal Tissue'))
result_TumorVsParatumor_OR <- results(dds_out, contrast=c('condition',
                                                          'Primary Tumor',
                                                          'Solid Tissue Normal'))
vsd <- vst(dds_out, blind=F)
a <- plotPCA(vsd, 'condition')
pdf('2020young\\TPN_PCA.pdf')
a
dev.off()

# DEGene screen.
DF_result_PN <- data.frame(GeneSymbol=result_ParatumorVsNormal_OR@rownames,
                           result_ParatumorVsNormal_OR@listData, stringsAsFactors=F)
DF_result_TP <- data.frame(GeneSymbol=result_TumorVsParatumor_OR@rownames,
                           result_TumorVsParatumor_OR@listData, stringsAsFactors=F)
TPN_cochange <- merge(DF_result_PN, DF_result_TP, by='GeneSymbol')
TPN_cochange_logFC_Pval <- data.frame(GeneSymbol=TPN_cochange$GeneSymbol,
                                      logFC_PN=TPN_cochange$log2FoldChange.x,
                                      logFC_TP=TPN_cochange$log2FoldChange.y,
                                      Pval_PN=TPN_cochange$padj.x,
                                      Pval_TP=TPN_cochange$padj.y,
                                      stringsAsFactors=F)
TPN_cochange_logFC_Pval <- TPN_cochange_logFC_Pval[(TPN_cochange_logFC_Pval$logFC_PN>0&
                                                      TPN_cochange_logFC_Pval$logFC_TP>0)|
                                                     (TPN_cochange_logFC_Pval$logFC_PN<0&
                                                        TPN_cochange_logFC_Pval$logFC_TP<0),]
TPN_cochange_logFC_Pval <- na.omit(TPN_cochange_logFC_Pval)
TPN_cochange_logFC_Pval$logFC.mean <- apply(TPN_cochange_logFC_Pval[,2:3], 1, mean)
TPN_cochange_logFC_Pval$padj.mean <- apply(TPN_cochange_logFC_Pval[,4:5],1, mean)
TPN_cochange_logFC_Pval <- TPN_cochange_logFC_Pval[order(TPN_cochange_logFC_Pval$logFC.mean,
                                                         decreasing=T),]

#Draw volcano plot
resultProcessing <- function(result, sigLevel, logFC){
  log2_padj <- data.frame(GeneSymbol=result$GeneSymbol,
                          log2FC=result$logFC.mean,
                          padj=result$padj.mean)
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
       xlab='log2FC', ylab='-log10(padj)', type='p',
       col=log2_padj$sig, pch=20, cex=0.8)
  #text(log2_padj[which(log2_padj$sig!=1),2],
  #     -log10(log2_padj[which(log2_padj$sig!=1),3]),
  #     log2_padj[which(log2_padj$sig!=1),1],
  #     col=log2_padj[which(log2_padj$sig!=1),4],
  #     pos=3, cex=0.6)
  abline(v=c(-logFCthreshold,logFCthreshold), h=-log10(sigThreshold), col='grey', lty=2)
  dev.off()
}
a <- resultProcessing(TPN_cochange_logFC_Pval, 0.01, 2)
volcanoPlot('2020young\\TPN_DEG_volcano1.pdf', a, 0.01, 2)
# Draw volcano plot end

# GSEA analysis.
library(clusterProfiler)
geneList <- TPN_cochange_logFC_Pval$logFC.mean
names(geneList) <- TPN_cochange_logFC_Pval$GeneSymbol
gmtfile <- read.gmt('D:\\Program Files\\GSEA_4.0.2\\gmtfile\\msigdb.v7.0.symbols.gmt')
GSEAresults <- GSEA(geneList, TERM2GENE=gmtfile)
GSEAresults_df <- GSEAresults@result
write.csv(GSEAresults_df, file='GSEAresults_df.csv')
GSEA_enrichID <- GSEAresults@result$ID
library(enrichplot)
pdf('2020young\\GSEAplots.pdf')
for (i in GSEA_enrichID) {
  A <- gseaplot2(GSEAresults, i, title=i)
  print(A)
}
dev.off()

# Significantly changed genes analysis.
TPN_cochange_logFC_Pval_filter <-
  TPN_cochange_logFC_Pval[(TPN_cochange_logFC_Pval$logFC_PN>1&
                             TPN_cochange_logFC_Pval$logFC_TP>1&
                             TPN_cochange_logFC_Pval$Pval_PN<0.05&
                             TPN_cochange_logFC_Pval$Pval_TP<0.05)|
                            (TPN_cochange_logFC_Pval$logFC_PN<(-1)&
                               TPN_cochange_logFC_Pval$logFC_TP<(-1)&
                               TPN_cochange_logFC_Pval$Pval_PN<0.05&
                               TPN_cochange_logFC_Pval$Pval_TP<0.05),]

# DE gene annotiation.
gene_info <- read.csv('D:/PublicData/Data_from_NCBI/Homo_sapiens.gene_info.csv',
                      stringsAsFactors=F, header=T)
gene_info <- data.frame(GeneSymbol=gene_info$Symbol,
                        description=gene_info$description,
                        Full_name_from_nomenclature_authority=gene_info$Full_name_from_nomenclature_authority)
TPN_cochange_logFC_Pval_filter_anno <- merge(TPN_cochange_logFC_Pval_filter,
                                             gene_info, by='GeneSymbol',
                                             all.x=T, sort=F)
write.csv(TPN_cochange_logFC_Pval_filter_anno, 
          file='TPN_cochange_logFC_Pval_filter_anno.csv')

# DE gene enrichment analysis.
library(clusterProfiler)
library(org.Hs.eg.db)
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

TPN_cochange_ensemblID <- merge(TPN_cochange_logFC_Pval_filter,
                                ensemblID_gene, by='GeneSymbol')
TPN_cochange_ensemblID <- TPN_cochange_ensemblID$EnsemblID
TPN_GO_KEGG_enrichment <-  GO_KEGG_enrichment(TPN_cochange_ensemblID,'TPN_cochange', 0.05, 0.05)

pdf(file='2020young\\TPN_GO_enrichment.pdf')
dotplot(TPN_GO_KEGG_enrichment$CC, title='TPN_CC_enrichment', showCategory=15)
dotplot(TPN_GO_KEGG_enrichment$BP, title='TPN_BP_enrichment', showCategory=15)
dotplot(TPN_GO_KEGG_enrichment$MF, title='TPN_MF_enrichment', showCategory=15)
dotplot(TPN_GO_KEGG_enrichment$KK, title='TPN_KK_enrichment', showCategory=15)
dev.off()

# WGCNA analysis for directional changed genes. ####################################################
# restart R and library(WGCNA), make WGCNA to be the first package in environment.
library(WGCNA)
enableWGCNAThreads(nThreads=12)
# Prepare the WGCNA analysis needed data.
library(DESeq2)
dataExpr <- assay(vsd)
dataExpr <- dataExpr[row.names(dataExpr)%in%TPN_cochange_logFC_Pval$GeneSymbol,]
dataExpr <- t(dataExpr)
dataTrait <- TcgaGtex_phenotype_with_3groups_outlierRemoved
library(reshape2)
dataTrait <- dcast(dataTrait, sample~condition)
dataTrait <- data.frame(sample=dataTrait$sample,
                        ifelse(is.na(dataTrait[,-1]), 0, 1),
                        stringsAsFactors=F)
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector = powers, verbose = 5)
pdf('2020young\\powerEstimate.pdf')
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
net <- blockwiseModules(dataExpr, power = sft$powerEstimate,
                        maxBlockSize = 30000,
                        TOMType = 'unsigned', minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs=TRUE, corType = 'pearson',
                        maxPOutliers=1, loadTOMs=TRUE,
                        saveTOMFileBase = 'tmp.tom',
                        verbose = 3)
table(net$colors)
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
pdf('2020young\\cluster_dendrogram.pdf')
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
MEs_col <- net$MEs
library(stringr)
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(net$MEs),"ME",""))))
MEs_col$sample <- row.names(MEs_col)
MEs_colpheno <- merge(MEs_col, dataTrait, by='sample')
row.names(MEs_colpheno) <- MEs_colpheno$sample
MEs_colpheno <- MEs_colpheno[,-1]
MEs_colpheno <- orderMEs(MEs_colpheno)
pdf('2020young\\ENplot.pdf')
plotEigengeneNetworks(MEs_colpheno, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,3),
                      marHeatmap = c(8,9,0,2), plotDendrograms = T, 
                      xLabelsAngle = 45)
dev.off()

load(net$TOMFiles[1], verbose=T)
TOM <- as.matrix(TOM)
dissTOM <- 1-TOM
plotTOM <- dissTOM^sft$powerEstimate
diag(plotTOM) <- NA
tomplot <- TOMplot(plotTOM, net$dendrograms, moduleColors, 
                   main = "Network heatmap plot, all genes")
save.image()

probes <- colnames(dataExpr)
dimnames(TOM) <- list(probes, probes)
cyt <- exportNetworkToCytoscape(TOM, edgeFile='cyt.edges.txt',
                                nodeFile='cyt.nodes.txt',
                                weighted=TRUE, threshold=0,
                                nodeNames=probes, nodeAttr=moduleColors)

dataTrait <- dataTrait[match(row.names(dataExpr), dataTrait$sample), ]
modTraitCor <- cor(MEs_col, dataTrait, use = "p")
modTraitP <- corPvalueStudent(modTraitCor, nrow(dataExpr))
textMatrix <- paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
modTraitCor <- modTraitCor[-6,-1]
textMatrix <- textMatrix[-6,-1]
textMatrix[5,2] <- '-0.92\n(1e-259)'
pdf('2020young\\labeledheatmap.pdf')
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(modTraitCor), 
               yLabels = row.names(modTraitCor), 
               cex.lab = 1, 
               ySymbols = row.names(modTraitCor), colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = TRUE, 
               cex.text = 1, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

geneModuleMembership = as.data.frame(cor(dataExpr, MEs_col, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(
  as.matrix(geneModuleMembership), ncol(dataExpr)))
geneTraitCor = as.data.frame(cor(dataExpr, dataTrait, use = "p"))
geneTraitP = as.data.frame(corPvalueStudent(
  as.matrix(geneTraitCor), ncol(dataExpr)))

module = "grey"
pheno = "Normal.Tissue"
modNames = substring(colnames(MEs_col), 3)
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(dataTrait))
moduleGenes = moduleColors == module
sizeGrWindow(7, 7)
pdf('2020young\\WGCNA.pdf')
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

save.image()


# WGCNA analysis for undirectional changed genes. ####################################################
# restart R and library(WGCNA), make WGCNA to be the first package in environment.
library(WGCNA)
enableWGCNAThreads(nThreads=12)
# Prepare the WGCNA analysis needed data.
library(DESeq2)
dataExpr <- assay(vsd)
dataExpr <- dataExpr[!row.names(dataExpr)%in%TPN_cochange_logFC_Pval$GeneSymbol,]
dataExpr <- t(dataExpr)
dataTrait <- TcgaGtex_phenotype_with_3groups_outlierRemoved
library(reshape2)
dataTrait <- dcast(dataTrait, sample~condition)
dataTrait <- data.frame(sample=dataTrait$sample,
                        ifelse(is.na(dataTrait[,-1]), 0, 1),
                        stringsAsFactors=F)
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector = powers, verbose = 5)
pdf('2020young\\powerEstimate_un.pdf')
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
net <- blockwiseModules(dataExpr, power = sft$powerEstimate,
                        maxBlockSize = ncol(dataExpr),
                        TOMType = 'unsigned', minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs=TRUE, corType = 'pearson',
                        maxPOutliers=1, loadTOMs=TRUE,
                        saveTOMFileBase = 'tmp.tom',
                        verbose = 3)
table(net$colors)
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
pdf('2020young\\cluster_dendrogram_un.pdf')
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
MEs_col <- net$MEs
library(stringr)
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(net$MEs),"ME",""))))
MEs_col$sample <- row.names(MEs_col)
MEs_colpheno <- merge(MEs_col, dataTrait, by='sample')
row.names(MEs_colpheno) <- MEs_colpheno$sample
MEs_colpheno <- MEs_colpheno[,-1]
MEs_colpheno <- orderMEs(MEs_colpheno)
pdf('2020young\\ENplot_un.pdf')
plotEigengeneNetworks(MEs_colpheno, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,3),
                      marHeatmap = c(8,9,0,2), plotDendrograms = T, 
                      xLabelsAngle = 45)
dev.off()

load(net$TOMFiles[1], verbose=T)
TOM <- as.matrix(TOM)
dissTOM <- 1-TOM
plotTOM <- dissTOM^sft$powerEstimate
diag(plotTOM) <- NA
tomplot <- TOMplot(plotTOM, net$dendrograms, moduleColors, 
                   main = "Network heatmap plot, all genes")

probes <- colnames(dataExpr)
dimnames(TOM) <- list(probes, probes)
cyt <- exportNetworkToCytoscape(TOM, edgeFile='cyt.edges.txt',
                                nodeFile='cyt.nodes.txt',
                                weighted=TRUE, threshold=0,
                                nodeNames=probes, nodeAttr=moduleColors)

dataTrait <- dataTrait[match(row.names(dataExpr), dataTrait$sample), ]
modTraitCor <- cor(MEs_col, dataTrait, use = "p")
modTraitP <- corPvalueStudent(modTraitCor, nrow(dataExpr))
textMatrix <- paste(signif(modTraitCor, 2), " (", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
modTraitCor <- modTraitCor[-dim(modTraitCor)[1],-1]
textMatrix <- textMatrix[-dim(textMatrix)[1],-1]
pdf('2020young\\labeledheatmap_un.pdf')
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(modTraitCor), 
               yLabels = row.names(modTraitCor), 
               cex.lab = 0.6, 
               ySymbols = row.names(modTraitCor), colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = TRUE, 
               cex.text = 0.6, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

geneModuleMembership = as.data.frame(cor(dataExpr, MEs_col, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(
  as.matrix(geneModuleMembership), ncol(dataExpr)))
geneTraitCor = as.data.frame(cor(dataExpr, dataTrait, use = "p"))
geneTraitP = as.data.frame(corPvalueStudent(as.matrix(geneTraitCor), ncol(dataExpr)))

module = "grey"
pheno = "Normal.Tissue"
modNames = substring(colnames(MEs_col), 3)
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(dataTrait))
moduleGenes = moduleColors == module
sizeGrWindow(7, 7)
pdf('2020young\\WGCNA_un.pdf')
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

save.image()

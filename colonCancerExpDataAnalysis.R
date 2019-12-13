
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
rm(GDC_colonCancer_ExpData)
GDC_colonCancer_ExpData_gene <- data.frame(GeneSymbol=GDC_colonCancer_ExpData_gene$GeneSymbol,
                                           GDC_colonCancer_ExpData_gene[,c(-1,-length(GDC_colonCancer_ExpData_gene))],
                                           stringsAsFactors=F)
GDC_colonCancer_ExpData_gene_agr <- aggregate(GDC_colonCancer_ExpData_gene[,-1],
                                              by=list(GDC_colonCancer_ExpData_gene$GeneSymbol),
                                              FUN=sum)
sample_id <- data.frame(sample_id=names(GDC_colonCancer_ExpData_gene_agr)[-1], stringsAsFactors=F)
rm(GDC_colonCancer_ExpData_gene)
GDC_colonCancer_ExpData_gene_agr_matrix <- as.matrix(x=GDC_colonCancer_ExpData_gene_agr)
row.names(GDC_colonCancer_ExpData_gene_agr_matrix) <- GDC_colonCancer_ExpData_gene_agr_matrix[,1]
GDC_colonCancer_ExpData_gene_agr_matrix <- GDC_colonCancer_ExpData_gene_agr_matrix[,-1]
GDC_colonCancer_ExpData_gene_agr_matrix <- apply(GDC_colonCancer_ExpData_gene_agr_matrix,
                                                 2, as.numeric)
# Above step will eliminate rownames of matrix. SHIT!
row.names(GDC_colonCancer_ExpData_gene_agr_matrix) <- GDC_colonCancer_ExpData_gene_agr$Group.1
rm(GDC_colonCancer_ExpData_gene_agr)
GDC_colonCancer_ExpData_gene_agr_matrix_filtered <-
  GDC_colonCancer_ExpData_gene_agr_matrix[rowSums(GDC_colonCancer_ExpData_gene_agr_matrix)!=0,]
save(GDC_colonCancer_ExpData_gene_agr_matrix_filtered,
     file='GDC_colonCancer_ExpData_gene_agr_matrix.RData')

cli_info <- read.csv('D:\\PublicData\\DataFromXena\\GDC_colon_clinical_info.csv',header=T,
                     stringsAsFactors=F)
stage_batch <- data.frame(sample_id=cli_info[,1], stage=cli_info$tumor_stage.diagnoses,
                          batch=cli_info$batch_number, stringsAsFactors=F)
stage_batch$sample_id <- gsub(x=stage_batch$sample_id, pattern='-', replacement='\\.')
group <- merge(sample_id, stage_batch, by='sample_id', all.x=T, sort=F)
group <- group[group$stage!=''&group$stage!='not reported',]
GDC_colonCancer_ExpData_gene_agr_matrix_filtered <-
  GDC_colonCancer_ExpData_gene_agr_matrix_filtered[,which(colnames(GDC_colonCancer_ExpData_gene_agr_matrix_filtered)%in%group$sample_id)]

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
dds <- DESeqDataSetFromMatrix(countData=GDC_colonCancer_ExpData_gene_agr_matrix_filtered,
                              colData=group, design=~batch+stage)
dds_out <- DESeq(dds, parallel=T)
save(dds_out, file='GDC_colonCancer_stage_dds_out.RData')
vsd <- vst(dds_out, blind=F)
plotPCA(vsd, 'stage')
assay(vsd) <- removeBatchEffect(assay(vsd), vsd$batch)
plotPCA(vsd, 'batch')
plotPCA(vsd, 'stage')

Stage2VsStage1 <- results(object=dds_out, contrast=c('stage', 'stage ii', 'stage i'))
Stage3VsStage2 <- results(object=dds_out, contrast=c('stage', 'stage iii', 'stage ii'))
stage4VsStage3 <- results(object=dds_out, contrast=c('stage', 'stage iv', 'stage iii'))
S2VsS1_DF <- data.frame(GeneSymbol=Stage2VsStage1@rownames, Stage2VsStage1@listData, stringsAsFactors=F)
S3VsS2_DF <- data.frame(GeneSymbol=Stage3VsStage2@rownames, Stage3VsStage2@listData, stringsAsFactors=F)
S4VsS3_DF <- data.frame(GeneSymbol=stage4VsStage3@rownames, stage4VsStage3@listData, stringsAsFactors=F)
tmp <- merge(S2VsS1_DF, S3VsS2_DF, by='GeneSymbol')
stage_results <- merge(tmp, S4VsS3_DF, by='GeneSymbol')
stage_logFC_Pval <- data.frame(GeneSymbol=stage_results$GeneSymbol,
                          logFC2Vs1=stage_results$log2FoldChange.x,
                          logFC3Vs2=stage_results$log2FoldChange.y,
                          logFC4Vs3=stage_results$log2FoldChange,
                          Pval_2Vs1=stage_results$padj.x,
                          Pval_3Vs2=stage_results$padj.y,
                          Pval_4Vs3=stage_results$padj,
                          stringsAsFactors=F)
stage_logFC_Pval_filter <- stage_logFC_Pval[(stage_logFC_Pval$logFC2Vs1>0&
                                               stage_logFC_Pval$logFC3Vs2>0&
                                               stage_logFC_Pval$logFC4Vs3>0&
                                               (stage_logFC_Pval$Pval_2Vs1<0.05|
                                               stage_logFC_Pval$Pval_3Vs2<0.05|
                                               stage_logFC_Pval$Pval_4Vs3<0.05))|
                                              (stage_logFC_Pval$logFC2Vs1<0&
                                                 stage_logFC_Pval$logFC3Vs2<0&
                                                 stage_logFC_Pval$logFC4Vs3<0&
                                                 (stage_logFC_Pval$Pval_2Vs1<0.05|
                                                 stage_logFC_Pval$Pval_3Vs2<0.05|
                                                 stage_logFC_Pval$Pval_4Vs3<0.05)),]
stage_logFC_Pval_filter <- na.omit(stage_logFC_Pval_filter)
stage_logFC_Pval_filter$logFC.mean <- apply(stage_logFC_Pval_filter[,2:4],1,mean)
stage_logFC_Pval_filter$padj.mean <- apply(stage_logFC_Pval_filter[,5:7],1,mean)
stage_logFC_Pval_filter_order <-
  stage_logFC_Pval_filter[order(abs(stage_logFC_Pval_filter$logFC.mean), decreasing=T),]

stage_logFC_Pval_filter_order_anno <- merge(stage_logFC_Pval_filter_order, gene_info, 
                                            by.x='GeneSymbol', by.y='Symbol', all.x=T, sort=F)
write.csv(stage_logFC_Pval_filter_order_anno, file='stage_logFC_Pval_filter_order_anno.csv')


stage_cochange_ensemblID <- merge(stage_logFC_Pval_filter_order, ensemblID_gene_GDC, by='GeneSymbol')
stage_cochange_ensemblID <- stage_cochange_ensemblID$Ensembl_ID
stage_GO_KEGG_enrichment <-  GO_KEGG_enrichment(stage_cochange_ensemblID,'stage_cochange', 0.05, 0.05)

write.csv(stage_GO_KEGG_enrichment$BP@result, file='stage_GO_KEGG_enrichment_BP.csv')
write.csv(stage_GO_KEGG_enrichment$CC@result, file='stage_GO_KEGG_enrichment_CC.csv')


png(file='stage_CC_enrichment.png')
dotplot(stage_GO_KEGG_enrichment$CC, title='stage_CC_enrichment', showCategory=15)
dev.off()
png(file='stage_BP_enrichment.png')
dotplot(stage_GO_KEGG_enrichment$BP, title='stage_BP_enrichment', showCategory=15)
dev.off()
dotplot(stage_GO_KEGG_enrichment$MF, title='stage_MF_enrichment', showCategory=15)
dotplot(stage_GO_KEGG_enrichment$KK, title='stage_KK_enrichment', showCategory=15)
dev.off()


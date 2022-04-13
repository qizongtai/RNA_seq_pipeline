library("DESeq2")
library("ggplot2")
library("RColorBrewer")
suppressMessages(library("gplots"))
library("pheatmap")
library("amap")
library("reshape2")
library("ggrepel")
# suppressMessages(library(DESeq2))
library("readr")
setwd("/scratch/splab/zqi/rna_seq/1_202009_tgfb_jq1_DEseq")
#rm(list=ls())

###----------read in matrix-------------###
#---count data
countData = read.table("/scratch/splab/zqi/rnaseq/1_data/all.gene_counts.tsv",header = T,sep = "\t", quote="")
#prepare data for DESeq2
rnames = make.names(countData$external_gene_name, unique = TRUE)
countData[,c(1:7)] = NULL
countData[is.na(countData)] <- 0
#countData <- as.matrix(sapply(countData, as.integer))
countData <- as.matrix(sapply(countData, function(x){round(x, digits = 0)}))
rownames(countData) <- rnames
# #reorder columns
# countData <- countData[,c("sample.scc9_jq1_r1",   "sample.scc9_jq1_r2",   "sample.scc9_jq1_r3")]
# #rename columns
colnames(countData) <- paste0(colnames(countData), c(rep("frozen", 11), rep("FFPE", 11))) %>% 
                       substr(.,8,length(.))
#

#---column data
colData = DataFrame(conditions = factor(c(rep("frozen", 11), rep("FFPE", 11))))
#---creat DESeq obj with count data and column data
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ conditions)
print(paste("Read in", nrow(dds),"genes"))
#3 replicates each for 12 conditions; nrow(colData)=total=3x12=36
keep <- rowSums(counts(dds)) > nrow(colData)/3
dds <- dds[keep,]
print(paste(nrow(dds),"genes remained after filtering of genes with all counts less than", nrow(colData)/2, "in all samples"))
dds <- DESeq(dds)
saveRDS(dds, file="3_output/rnaseq_36samples_count_morethan12_deseq2_dds.rds")

###--------read from RDS--------###
# dds <- readRDS(file = "3_output/rnaseq_36samples_count_morethan12_deseq2_dds.rds")
# str(dds)

###-------settings-------####
output_prefix <- "/scratch/splab/zqi/rnaseq/3_output"

###-------Get normalized counts-------###
print("Output normalized counts")
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
normalized_counts_output = data.frame(id=rownames(normalized_counts), normalized_counts)
write.table(normalized_counts_output, file=paste0(output_prefix,".DESeq2.normalized.xls"),
            quote=F, sep="\t", row.names=F, col.names=T)

norm_mat_melt <- reshape2::melt(normalized_counts_output, id.vars = c('id'))
pdf(file=paste0(output_prefix,".DESeq2.normalized.norm.pdf"), width=12, height=7, useDingbats=FALSE)
ggplot(norm_mat_melt, aes(x=variable, y=value)) + geom_boxplot(aes(color=variable)) +
  geom_violin(aes(fill=variable), alpha=0.5) + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank()) + ylab("norm transformed expression value")
dev.off()

###-------Get transformed counts-------###
print("Output transformed counts")
#rlog transform
rld <- rlog(dds, blind=FALSE)
rlogMat <- assay(rld)
rlogMat <- rlogMat[order(normalized_counts_mad, decreasing=T), ]
print("Output rlog transformed normalized counts")
rlogMat_output = data.frame(id=rownames(rlogMat), rlogMat)
write.table(rlogMat_output, file=paste0(output_prefix,".DESeq2.normalized.rlog.xls"),
            quote=F, sep="\t", row.names=F, col.names=T)

rlog_mat_melt <- reshape2::melt(rlogMat_output, id.vars = c('id'))
pdf(file=paste0(output_prefix,".DESeq2.normalized.rlog.pdf"), width=12, height=7, useDingbats=FALSE)
ggplot(rlog_mat_melt, aes(x=variable, y=value)) + geom_boxplot(aes(color=variable)) +
  geom_violin(aes(fill=variable), alpha=0.5) + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank()) + ylab("rLog transformed expression value")
dev.off()

#vst transform
vsd <- vst(dds, blind=FALSE)
vstMat <- assay(vsd)
vstMat <- vstMat[order(normalized_counts_mad, decreasing=T), ]
print("Output vst transformed normalized counts")
vstMat_output = data.frame(id=rownames(vstMat), vstMat)
write.table(vstMat_output, file=paste0(output_prefix,".DESeq2.normalized.vst.xls"),
            quote=F, sep="\t", row.names=F, col.names=T)

vst_mat_melt <- reshape2::melt(vstMat_output, id.vars = c('id'))
pdf(file=paste0(output_prefix,".DESeq2.normalized.vst.pdf"), width=12, height=7, useDingbats=FALSE)
ggplot(vst_mat_melt, aes(x=variable, y=value)) + geom_boxplot(aes(color=variable)) +
  geom_violin(aes(fill=variable), alpha=0.5) + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank()) + ylab("vst transformed expression value")
dev.off()

###-------sample level analysis: hc clustering-------###
print("Performing sample clustering")
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
pearson_cor <- round(as.matrix(cor(rlogMat, method="pearson")),4)

hc <- hcluster(t(rlogMat), method="pearson")

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

pearson_cor <- pearson_cor[hc$order, hc$order]

pearson_cor_output = data.frame(id=rownames(pearson_cor), pearson_cor)
write.table(pearson_cor_output, file=paste0(output_prefix,".DESeq2.pearson_cor.xls"),
            quote=F, sep="\t", row.names=F, col.names=T)

upper_tri <- get_upper_tri(pearson_cor)
# Melt the correlation matrix
melted_cormat <- reshape2::melt(upper_tri, na.rm = TRUE)

col = colorRampPalette(rev(brewer.pal(n=7, name="RdYlBu")))(100)

# Create a ggheatmap
p <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradientn(colours=col, name="Pearson correlation") + theme_classic() +
  coord_fixed() + 
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.justification = c(1, 0),
    legend.position = "top",
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 8, barheight = 1, title.position = "left"))

pdf(file=paste0(output_prefix,".DESeq2.normalized.rlog.pearson.pdf"), width=9, height=9, useDingbats=FALSE)
p
dev.off()
#ggsave(filename=paste0(output_prefix,".DESeq2.normalized.rlog.pearson.pdf"),width=13.5,height=15,units=c("cm"))

###-------sample level analysis: PCA-------###
#method1: baseR fun prcomp
print("PCA analysis")
formulaV <- c("conditions")
topn = 5000
rlogMat_nrow = nrow(rlogMat)
if (topn > rlogMat_nrow){
  topn = rlogMat_nrow
}
pca_mat = rlogMat[1:topn,]
pca_mat <- as.data.frame(t(pca_mat))
pca <- prcomp(pca_mat, scale=T)
pca_x = pca$x
pca_individual = data.frame(samp=rownames(pca_x), pca_x, colData)
write.table(pca_individual, file=paste0(output_prefix,".DESeq2.pca_individuals.xls"), sep="\t", quote=F, row.names=F, col.names=T)
pca_percentvar <- formatC(pca$sdev^2 * 100 / sum( pca$sdev^2))

if (length(formulaV)==1) {
  p <- ggplot(pca_individual, aes(PC1, PC2, color=conditions))
} else if (length(formulaV==2)) {
  p <- ggplot(pca_data, aes(PC1, PC2, color=conditions,shape=conditions))
}
p = p + geom_point(size=3) + 
  xlab(paste0("PC1: ", pca_percentvar[1], "% variance")) +
  ylab(paste0("PC2: ", pca_percentvar[2], "% variance")) +
  geom_text_repel(aes(label=samp), show.legend=F) +
  theme_classic() +
  theme(legend.position="top", legend.title=element_blank())

pdf(file=paste0(output_prefix,".DESeq2.normalized.rlog.pca.pdf"), width=10, height=10, useDingbats=FALSE)
p
dev.off()
#ggsave(p, filename=paste0(output_prefix,".DESeq2.normalized.rlog.pca.pdf"),width=10,height=10,units=c("cm"))

pca_percentvar <- data.frame(PC=colnames(pca_x), Variance=pca_percentvar)
write.table(pca_percentvar, file=paste0(output_prefix,".DESeq2.pca_pc_weights.xls"), sep="\t", quote=F, row.names=F, col.names=T)

#method2: DESeq2::plotPCA
plotPCA(rld, intgroup = c("conditions"))


###-------gene level analysis-------###
###-------two group comparison------###
if (file.exists(paste0(output_prefix,".DESeq2.all.DE"))) 
file.remove(paste0(output_prefix,".DESeq2.all.DE"))
# sampleA <- "scc47_tgfb"
# sampleB <- "scc47_con"

# sampleA <- "scc4_tgfb"
# sampleB <- "scc4_con"

# sampleA <- "scc9_tgfb"
# sampleB <- "scc9_con"

# sampleA <- "jhu6_tgfb"
# sampleB <- "jhu6_con"

#-----------------------
sampleA <- "scc47_jq1"
sampleB <- "scc47_con"

# sampleA <- "scc4_jq1"
# sampleB <- "scc4_con"

# sampleA <- "scc9_jq1"
# sampleB <- "scc9_con"

# sampleA <- "jhu6_jq1"
# sampleB <- "jhu6_con"

# no changes are needed for the below code
print(paste("DE genes between", sampleA, sampleB, sep=" "))
contrastV <- c("conditions", sampleA, sampleB)
res <- results(dds,  contrast=contrastV)
baseA <- counts(dds, normalized=TRUE)[, colData(dds)$conditions == sampleA]
if (is.vector(baseA)){
  baseMeanA <- as.data.frame(baseA)
} else {
  baseMeanA <- as.data.frame(rowMeans(baseA))
}
baseMeanA <- round(baseMeanA, 3)
colnames(baseMeanA) <- sampleA
baseB <- counts(dds, normalized=TRUE)[, colData(dds)$conditions == sampleB]
if (is.vector(baseB)){
  baseMeanB <- as.data.frame(baseB)
} else {
  baseMeanB <- as.data.frame(rowMeans(baseB))
}
baseMeanB <- round(baseMeanB, 3)
colnames(baseMeanB) <- sampleB
res <- cbind(baseMeanA, baseMeanB, as.data.frame(res))
res <- data.frame(ID=rownames(res), res)
res$baseMean <- round(rowMeans(cbind(baseA, baseB)),3)
res$padj[is.na(res$padj)] <- 1
res$pvalue[is.na(res$pvalue)] <- 1
res$log2FoldChange <- round(res$log2FoldChange,3)
res$padj <- as.numeric(formatC(res$padj))
res$pvalue <- as.numeric(formatC(res$pvalue))

res <- res[order(res$padj),]
comp314 <- paste(sampleA, "_vs_", sampleB, sep=".")
file_base <- paste(output_prefix, "DESeq2", comp314, sep=".")
file_base1 <- paste(file_base, "results.xls", sep=".")
res_output <- as.data.frame(subset(res,select=c('ID',sampleA,sampleB,"baseMean",'log2FoldChange','pvalue', 'padj')))
write.table(res_output, file=file_base1, sep="\t", quote=F, row.names=F)

res_de <- subset(res, res$padj<0.05, select=c('ID', sampleA, sampleB, 'log2FoldChange', 'padj'))
res_de_up <- subset(res_de, res_de$log2FoldChange>=1)
file <- paste(output_prefix, "DESeq2",sampleA, "_higherThan_", sampleB, 'xls', sep=".") 
write.table(as.data.frame(res_de_up), file=file, sep="\t", quote=F, row.names=F)
res_de_up_id <- subset(res_de_up, select=c("ID"))
#file <- paste(file_base, "DE_up_id", sep=".")
file <- paste(output_prefix, "DESeq2", sampleA, "_higherThan_", sampleB,'id.xls', sep=".") 
write.table(as.data.frame(res_de_up_id), file=file, sep="\t", quote=F, row.names=F, col.names=F)

if(dim(res_de_up_id)[1]>0) {
  res_de_up_id_l <- cbind(res_de_up_id, paste(sampleA, "_higherThan_",sampleB, sep="."))
  write.table(as.data.frame(res_de_up_id_l), file=paste0(output_prefix,".DESeq2.all.DE"),
              sep="\t",quote=F, row.names=F, col.names=F, append=T)
}

res_de_dw <- subset(res_de, res_de$log2FoldChange<=(-1)*1)
file <- paste(output_prefix, "DESeq2",sampleA, "_lowerThan_", sampleB, 'xls', sep=".") 
write.table(as.data.frame(res_de_dw), file=file, sep="\t", quote=F, row.names=F)
res_de_dw_id <- subset(res_de_dw, select=c("ID"))
file <- paste(output_prefix, "DESeq2",sampleA, "_lowerThan_", sampleB, 'id.xls', sep=".") 
write.table(as.data.frame(res_de_dw_id), file=file, sep="\t", quote=F, row.names=F, col.names=F)

if(dim(res_de_dw_id)[1]>0) {
  res_de_dw_id_l <- cbind(res_de_dw_id, paste(sampleA, "_lowerThan_",sampleB, sep="."))
  write.table(as.data.frame(res_de_dw_id_l), file=paste0(output_prefix,".DESeq2.all.DE"),
              sep="\t",quote=F, row.names=F, col.names=F, append=T)
}

res_output$level <- ifelse(res_output$padj<=0.05, ifelse(res_output$log2FoldChange>=1, paste(sampleA,"UP"), ifelse(res_output$log2FoldChange<=1*(-1), paste(sampleB,"UP"), "NoDiff")) , "NoDiff")
res_output$padj <- (-1) * log10(res_output$padj)
res_output$padj <- replace(res_output$pad, res_output$pad>5, 5.005)

boundary = ceiling(max(abs(res_output$log2FoldChange)))
p = ggplot(res_output, aes(x=log2FoldChange,y=padj,color=level)) +
  #geom_point(aes(size=baseMean), alpha=0.5) + theme_classic() +
  geom_point(size=1, alpha=0.5) + theme_classic() +
  xlab("Log2 transformed fold change") + ylab("Negative Log10 transformed FDR") +
  xlim(-1 * boundary, boundary) + theme(legend.position="top", legend.title=element_blank())

pdf(file=paste0(file_base1,".volcano.pdf"), width=8, height=8, useDingbats=FALSE)
p
dev.off()
#ggsave(p, filename=paste0(file_base1,".volcano.pdf"),width=13.5,height=15,units=c("cm"))

#------visualization1: gene ranking by logFC------#
#head(res_output)
#lable gene names that have log2FC > 5
res_output_line <- res_output[order(res_output$log2FoldChange),]
res_output_line$x <- 1:nrow(res_output_line)
res_output_line[abs(res_output_line$log2FoldChange)<5 | res_output_line$level=="NoDiff", "ID"] <- NA
head(res_output_line)

res_output_line[complete.cases(res_output_line), ]

library(ggrepel)
p <- ggplot(res_output_line, aes(x=x, y=log2FoldChange)) + geom_point(aes(color=log2FoldChange)) + 
  scale_color_gradient2(low="dark green", mid="yellow", high= "dark red", midpoint = 0) + 
  theme_classic() + geom_hline(yintercept = 0, linetype="dotted")
# label gene names
p <- p +  geom_text_repel(aes(label=ID))
pdf(file=paste0(output_prefix, ".", sampleA, sampleB, ".de.log2fc.morethan5.pdf"), width=8, height=8, useDingbats=FALSE)
p
dev.off()

#------visualization2: heatmap------#
res_de_up_top20_id <- as.vector(head(res_de_up$ID,20))
res_de_dw_top20_id <- as.vector(head(res_de_dw$ID,20))

res_de_top20 <- c(res_de_up_top20_id, res_de_dw_top20_id)
res_de_top20

res_de_top20_expr <- normalized_counts[res_de_top20,]
#library(pheatmap)
#pheatmap(res_de_top20_expr, cluster_row=T, scale="row", annotation_col=sample)

#for scc47 tgfb vs con
res_de_top20_expr_subset <- res_de_top20_expr[,c(13:18)]
# #for scc4 tgfb vs con
# res_de_top20_expr_subset <- res_de_top20_expr[,c(31:36)]
# #for scc9 tgfb vs con
# res_de_top20_expr_subset <- res_de_top20_expr[,c(4:9)]
# #for jhu6 tgfb vs con
# res_de_top20_expr_subset <- res_de_top20_expr[,c(22:27)]

library(ComplexHeatmap)
#p <- Heatmap(res_de_top20_expr_subset, cluster_rows = FALSE)
p <- Heatmap(log10(res_de_top20_expr_subset+0.1), cluster_rows = FALSE)
pdf(file=paste0(output_prefix, ".", sampleA, sampleB, ".top20.up.and.down.pdf"), width=8, height=12, useDingbats=FALSE)
p
dev.off()

#------visualization3: volcano plot------#
# add gene name column to df
res_de_top20_expr2 <- data.frame(Gene=rownames(res_de_top20_expr), res_de_top20_expr)
head(res_de_top20_expr2)
# to long format
res_de_top20_expr2 <- melt(res_de_top20_expr, id=c("Gene"))
colnames(res_de_top20_expr2) <- c("Gene", "Sample", "Expression")
head(res_de_top20_expr2)
# volcano plot
ggplot(res_de_top20_expr2, aes(x=Gene, y=Expression)) + 
  geom_point(aes(color=Sample), alpha=0.5) + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank()) + ylab("Normalized xpression value") + scale_y_log10()

# try with res_de_top20_expr_subset
# add gene name column to df
res_de_top20_expr2 <- data.frame(Gene=rownames(res_de_top20_expr_subset), res_de_top20_expr_subset)
head(res_de_top20_expr2)
# to long format
res_de_top20_expr2 <- melt(res_de_top20_expr_subset, id=c("Gene"))
colnames(res_de_top20_expr2) <- c("Gene", "Sample", "Expression")
head(res_de_top20_expr2)
# volcano plot
ggplot(res_de_top20_expr2, aes(x=Gene, y=Expression)) + 
  geom_point(aes(color=Sample), alpha=0.5) + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank()) + ylab("Normalized xpression value") + scale_y_log10()

###----prepare logFC ranked gene for GSEA------###
#---for tgfb genes
suppressMessages(library(data.table))
cell <- c("scc47", "jhu6", "scc4", "scc9")
for (i in cell) {
  input_suffix <- paste0(".DESeq2.", i ,"_tgfb._vs_.", i, "_con.results.xls")
  output_suffix <- paste0(".DESeq2.", i ,"_tgfb._vs_.", i, "_con.log2fc_ranked.symbol.rnk")
  vs_result <- read.table(paste0(output_prefix, input_suffix), header=T, sep="\t", quote="")
  vs_result <- subset(vs_result, select=c("ID","log2FoldChange"))
  
  # ###---this is to replace ensumbl gene ID with gene symble 
  # vs_result <- data.table(vs_result, key="ID")
  # idmap <- read.table("13_salmon_deseq2/genome/GRCh38.idmap", header=T, sep="\t", quote="")
  # idmap <- data.table(idmap, key="ENSG")
  # merge_result <- merge(vs_result, idmap, by.x="ID", by.y="ENSG", all.x=T)
  # merge_result <- merge_result[order(merge_result$log2FoldChange, decreasing = T),]
  # merge_result <- merge_result[merge_result$Symbol!="", c(3,2)]
  # write.table(merge_result, "13_salmon_deseq2/ehbio_salmon.DESeq2.log2fc_ranked.symbol.rnk", col.names = T, row.names = F, sep="\t", quote=F)
  
  vs_result <- vs_result[order(vs_result$log2FoldChange, decreasing = T),]
  vs_result <- vs_result[vs_result$ID!="", ]
  write.table(vs_result, paste0(output_prefix, output_suffix), col.names = T, row.names = F, sep="\t", quote=F)
}
#---for jq1 genes
suppressMessages(library(data.table))
cell <- c("scc47", "jhu6", "scc4", "scc9")
for (i in cell) {
  input_suffix <- paste0(".DESeq2.", i ,"_jq1._vs_.", i, "_con.results.xls")
  output_suffix <- paste0(".DESeq2.", i ,"_jq1._vs_.", i, "_con.log2fc_ranked.symbol.rnk")
  vs_result <- read.table(paste0(output_prefix, input_suffix), header=T, sep="\t", quote="")
  vs_result <- subset(vs_result, select=c("ID","log2FoldChange"))
  
  # ###---this is to replace ensumbl gene ID with gene symble 
  # vs_result <- data.table(vs_result, key="ID")
  # idmap <- read.table("13_salmon_deseq2/genome/GRCh38.idmap", header=T, sep="\t", quote="")
  # idmap <- data.table(idmap, key="ENSG")
  # merge_result <- merge(vs_result, idmap, by.x="ID", by.y="ENSG", all.x=T)
  # merge_result <- merge_result[order(merge_result$log2FoldChange, decreasing = T),]
  # merge_result <- merge_result[merge_result$Symbol!="", c(3,2)]
  # write.table(merge_result, "13_salmon_deseq2/ehbio_salmon.DESeq2.log2fc_ranked.symbol.rnk", col.names = T, row.names = F, sep="\t", quote=F)
  
  vs_result <- vs_result[order(vs_result$log2FoldChange, decreasing = T),]
  vs_result <- vs_result[vs_result$ID!="", ]
  write.table(vs_result, paste0(output_prefix, output_suffix), col.names = T, row.names = F, sep="\t", quote=F)
}

#############################---END---#######################################
# Filtering to remove rows with 0 reads
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)
res <- results(dds,contrast=c("condition","treated/CPI","untreated/DMSO"))
summary(res)
plotMA(res, main = "DESeq2,Female CPI versus DMSO")
resOrdered <- res[order(res$padj),]
write.table(resOrdered,"mitra_lab/CPI_RNAseq/3_output/f_all_DMSOvsCPI.txt",quote = F,sep ="\t")
sig <- resOrdered[!is.na(resOrdered$padj) &
                    resOrdered$padj<0.10 &
                    abs(resOrdered$log2FoldChange) >= fold,]
write.table(sig,"mitra_lab/CPI_RNAseq/3_output/f_sig_DMSOvsCPI.txt",quote = F,sep ="\t")
sig <- resOrdered[!is.na(resOrdered$padj) &
                    resOrdered$padj<0.10 &
                    resOrdered$log2FoldChange >= fold,]
write.table(sig,"mitra_lab/CPI_RNAseq/3_output/f_sigdown_DMSOvsCPI.txt",quote = F,sep ="\t")
sig <- resOrdered[!is.na(resOrdered$padj) &
                    resOrdered$padj<0.10 &
                    resOrdered$log2FoldChange <= -fold,]
write.table(sig,"mitra_lab/CPI_RNAseq/3_output/f_sigup_DMSOvsCPI.txt",quote = F,sep ="\t")
sig <- resOrdered[!is.na(resOrdered$padj) &
                    resOrdered$padj>=0.10,]
write.table(sig,"mitra_lab/CPI_RNAseq/3_output/f_notsig_DMS.txt",quote = F,sep ="\t")

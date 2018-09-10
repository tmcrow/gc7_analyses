source("https://bioconductor.org/biocLite.R")
biocLite("voom")
library("edgeR")
sample_submission <- read.csv("~/Box Sync/inv_4_growthchamber/libraries/sample_submission/sample_submission_new.csv")
sample_submission$inv4m_genotype <- as.character(sample_submission$inv4m_genotype)
sample_submission <- sample_submission[sample_submission$line %in% c("PT_NIL", "Mi21_NIL"),]
sample_submission <- sample_submission[sample_submission$inv4m_genotype %in% c("B73", "Inv4m"),]
sample_submission <- sample_submission[,c("tissue","temp_tx","inv4m_genotype","samp")]
sample_submission$samp <- as.character(sample_submission$samp)

gene_cnts <- read.csv("~/Box Sync/Taylor/gc7_data/gc7_tx_counts.csv")
gc <- gene_cnts[,names(gene_cnts) %in% sample_submission$samp]

sample_submission <- sample_submission[which(sample_submission$samp %in% names(gc)),]

cpm_log <- cpm(gc, log = TRUE)
median_log2_cpm <- apply(cpm_log, 1, median)
hist(median_log2_cpm)
expr_cutoff <- -1
abline(v = expr_cutoff, col = "red", lwd = 3)
sum(median_log2_cpm > expr_cutoff)
data_clean <- gc[median_log2_cpm > expr_cutoff, ]
cpm_log <- cpm(data_clean, log = TRUE)
heatmap(cor(cpm_log))

#
dim(t(cpm_log))
pca <- prcomp(t(cpm_log), scale. = TRUE)
tmp <- summary(pca)
length(unique(sample_submission$tissue))

require(RColorBrewer)
display.brewer.all()
str(sample_submission)
(jColors <- with(sample_submission,
         data.frame(tissue = levels(tissue),
                    color = I(brewer.pal(nlevels(tissue), name = 'Set1'))))) 
?display.brewer.pal
palette(brewer.pal(n = 9, name = 'Accent'))

plot(pca$x[, 1], pca$x[,2], pch=16, xlab=paste("PC1: ",tmp$importance[2]*100,"% ","variance explained", sep=""), 
     ylab=paste("PC2: ", tmp$importance[5]*100,"%"," variance explained", sep=""), col=sample_submission$tissue)
     
#
group <- sample_submission$inv4m_genotype
group
y <- DGEList(counts = data_clean, group = group)
y
#
y <- calcNormFactors(y)
y$samples
#
y <- estimateDisp(y)
sqrt(y$common.dispersion) # biological coefficient of variation
plotBCV(y)
#
cpm(gc)
head(gc)
norm_gc_tm <- calcNormFactors(gc, method="TMM")
norm_gc_rle <- calcNormFactors(gc, method="RLE")
norm_gc_uq <- calcNormFactors(gc, method="upperquartile")


# VOOM --------------------------------------------------------------------

require(voom)


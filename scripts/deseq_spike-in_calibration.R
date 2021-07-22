# Differential Gene Expression Analysis using DESeq
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Includes spike-in (dm6) normalisation and batch effect correction

## Output:
### Spike-in-based size factors
### PCA plot
### Heatmap of sample distances
### Table with normalised counts
### Table with calculated RPKM per sample
### Table with calculated RPKM per condition
### Differential results file
### Summary file (with 0.05 FDR cutoff; no LFC cutoff)
### MA plots

library(DESeq2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(vsn)
library(ggrastr)


# Defining parameters -----------------------------------------------------

## Name of project
p <- "PRC2degron.cnRNAseq"

## Conditions & replicates:
expConditions <- c("UNT", "dTAG")
expReplicates <- c("rep1", "rep2", "rep3")


# Preparing spike-in read counts table -----------------------------------------

## read counts from a set of unique dm6 refGene genes prepared using a SAMtools-based custom Perl script

dm6 <- read.delim("dm6.refGene_GeneBody_Uniq.ANNOTATED.txt", header = T)

names(dm6) <- gsub(pattern = "_ReadCount", replacement = "", x = names(dm6))
rownames(dm6) <- make.names(c(as.character(dm6$ID)), unique = TRUE)
dm6 <- dplyr::select(dm6, contains("rep"))

# Preparing spike-in sample information table -----------------------------

dm6.Condition <- rep(expConditions, times = 3)
dm6.Rep <- rep(expReplicates, each = 2)
dm6.SampleInfo <- data.frame(condition = dm6.Condition, rep = dm6.Rep)
rownames(dm6.SampleInfo) <- colnames(dm6)
dm6.SampleInfo
## explicitly setting the factors levels
dm6.SampleInfo$condition <- factor(dm6.SampleInfo$condition,
  levels = c("UNT", "dTAG")
)


# Obtaining size factors from dm6 read counts -----------------------------

## Measuring the effect of the condition, while controlling for batch (rep) effects:
dm6.DESEQ2 <- DESeqDataSetFromMatrix(
  countData = dm6,
  colData = dm6.SampleInfo,
  design = ~ rep + condition
)

cds_dm6 <- estimateSizeFactors(dm6.DESEQ2)
size_Factors_dm6 <- sizeFactors(cds_dm6)
write.table(size_Factors_dm6, paste0("SizeFactors_", p, ".txt"),
  quote = FALSE, col.names = F
)


# Preparing mm10 (experimental) read counts table --------------------------

## read counts from a non-redundant mm10 refGene gene set prepared using a SAMtools-based custom Perl script

exp <- read.delim("mm10.NonRedundantRefGene_GeneBody.ANNOTATED.txt", header = T)
names(exp) <- gsub(pattern = "_ReadCount", replacement = "", x = names(exp))

geneSize <- dplyr::select(exp, ID, Size) %>%
  rename(refID = ID)

rownames(exp) <- make.names(c(as.character(exp$ID)), unique = TRUE)
exp <- dplyr::select(exp, contains("rep"))


# Preparing mm10 sample info table ----------------------------------------

exp.Condition <- factor(rep(expConditions, times = 3))
exp.Rep <- factor(rep(expReplicates, each = 2))
exp.SampleInfo <- data.frame(condition = exp.Condition, rep = exp.Rep)
rownames(exp.SampleInfo) <- colnames(exp)
exp.SampleInfo
## explicitly setting the factors levels
exp.SampleInfo$condition <- factor(exp.SampleInfo$condition,
  levels = c("UNT", "dTAG")
)

# Performing Differential Expression Analysis -----------------------------

## With variable of interest (condition) at the end of the design formula
exp.DESEQ2 <- DESeqDataSetFromMatrix(
  countData = exp,
  colData = exp.SampleInfo,
  design = ~ rep + condition
)

## Setting size factors to spike-in ones
sizeFactors(exp.DESEQ2) <- size_Factors_dm6


## Differential expression analysis:
exp.DESEQ2.analysis <- DESeq(exp.DESEQ2)

# Count data transformation ----------------------------------

### using variance stabilizing transformation (VST) ->
### making count values approximately homoskedastic (having constant variance along the range of mean values)
### If many of genes have large differences in counts due to the experimental design, it is important to set blind=FALSE

vsd <- vst(exp.DESEQ2.analysis, blind = F)
meanSdPlot(assay(vsd))

## Calculating Euclidean distance for the heatmap of the sample distances

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)


## Plotting heatmap of distance matrix
### gives an overview over similarities and dissimilarities between samples

colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pheatmap(sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  col = colors,
  main = "Euclidean distance sample clustering (VST-transformed data)",
  filename = paste0(p, ".SampleDistanceHeatmap.VST.pdf")
)

dev.off()

## Plotting PCA

pcaData <- plotPCA(vsd, intgroup = c("condition", "rep"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = condition, shape = rep)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ylim(-10, 10) +
  theme_minimal() +
  ggtitle("VST-transformed")

ggsave(paste0(p, ".PCA.VST.pdf"))


# Extracting normalised counts --------------------------------------------

normCounts <- as.data.frame(counts(exp.DESEQ2, normalized = TRUE))
normCounts <- rownames_to_column(normCounts, var = "refID")

normCounts <- right_join(geneSize, normCounts, by = "refID")

write.csv(normCounts, paste0(p, ".DESeq2_normCounts.csv"),
  quote = FALSE, row.names = F
)



# Calculating RPKMs  -----------------------

## Individual samples

rpkm <- mutate_at(normCounts, vars(contains("rep")), list(~ . / (Size / 1000)))
write.csv(rpkm, paste0(p, ".DESeq2_individualRPKMs.csv"),
  quote = FALSE, row.names = F
)

## Conditions means

rpkm.mean <- mutate(rpkm,
  UNT_rpkm = rowMeans(dplyr::select(rpkm, contains("UNT"))),
  dTAG_rpkm = rowMeans(dplyr::select(rpkm, contains("dTAG")))
) %>%
  dplyr::select(-contains("rep"))

write.csv(rpkm.mean, paste0(p, ".DESeq2_meanRPKMs.csv"),
  quote = FALSE, row.names = F
)


# Extracting results -------------------------------


## Raw LFC values

res <- results(exp.DESEQ2.analysis,
  name = "condition_dTAG_vs_UNT",
  alpha = 0.05
)

## plot default MA

pdf(paste0(p, ".builtinMAplot.pdf"))
plotMA(res, alpha = 0.05, colSig = "red")
dev.off()


## Shrinking LFC by method "apeglm"

res.apeglm <- lfcShrink(exp.DESEQ2.analysis,
  coef = "condition_dTAG_vs_UNT",
  res = res,
  type = "apeglm"
)

res$LFC_apeglm <- res.apeglm$log2FoldChange

## Exporting results with raw and shrunken log2FoldChange (LFC)

res.df <- as.data.frame(res) %>%
  rownames_to_column(var = "refID")

res.df <- right_join(geneSize, res.df, by = "refID") %>% 
  left_join(select(rpkm.mean, refID, UNT_rpkm), by = "refID") %>% 
  arrange(desc(LFC_apeglm))

write.csv(res.df,
  paste0(p, ".DESeq2_results.2hdTAG_vs_UNT.csv"),
  quote = F, row.names = F
)

## Summary of changes (w/o any LFC cutoff!)

sink(file = paste0(p, ".DESeq2_summary.txt"))
summary(res, alpha = 0.05)
sink()


# Making custom MA plots ----------------------------------------------------------

my_theme <- theme_bw() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
    axis.text = element_text(size = 16, colour = "black"),
    axis.title = element_text(size = 18, colour = "black", margin = 4),
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.position = "none"
  )


## Defining UP and DOWN genes based on padj cutoff and shrunken LFC cutoff

res.df$Reg_apeglm <- "n.s."
res.df[which(res.df$padj < 0.05 & res.df$LFC_apeglm > log2(1.5)), ]$Reg_apeglm <- "UP"
res.df[which(res.df$padj < 0.05 & res.df$LFC_apeglm < -log2(1.5)), ]$Reg_apeglm <- "DOWN"
table(res.df$Reg_apeglm)


res.df$Reg_norm <- "n.s."
res.df[which(res.df$padj < 0.05 & res.df$LFC_norm > log2(1.5)), ]$Reg_norm <- "UP"
res.df[which(res.df$padj < 0.05 & res.df$LFC_norm < -log2(1.5)), ]$Reg_norm <- "DOWN"
table(res.df$Reg_norm)



# Plotting MA with apeglm-shrunken LFC values -----------------------------

## numbers for the plot title:
UP <- nrow(res.df[res.df$Reg_apeglm == "UP", ])
DOWN <- nrow(res.df[res.df$Reg_apeglm == "DOWN", ])

ggplot(res.df, aes(log2(UNT_rpkm), LFC_apeglm)) +
  geom_point_rast(colour = "darkgrey",alpha = 0.3, raster.dpi = 600, size = 2) +
  geom_point_rast(data = res.df[res.df$Reg_apeglm != "n.s.", ],
                  aes(log2(UNT_rpkm), LFC_apeglm),
                  colour = "red", alpha = 0.3, raster.dpi = 600, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1, colour = "black") +
  xlab("log2 UNT RPKM") +
  ylab("log2 Fold Change (2h dTAG / UNT)") +
  ggtitle(paste0("UP = ", UP, ", DOWN = ", DOWN)) +
  my_theme +
  scale_x_continuous(limits = c(-10, 18)) +
  scale_y_continuous(limits = c(-3, 3))

ggsave(paste0(p, ".DESeq2_MAplot.pdf"),
       width = 5.2, height = 6)
# For Figure 2d, making heatmaps depicting gene expression (cnRNA-seq) in untreated cells
# and changes following IAA treatment for the three groups of genes defined by
# the earliest time of derepression. Heatmaps are sorted by RING1B cChIP-seq signal as in Figure 2e.

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)


# Importing data ----------------------------------------------------------

## FirstTimeUp gene grouping defined by the earliest time of derepression, ordered by RING1B signal in 5 kb promoters (deepTools)
first.order <- read.delim("../data/FirstTimeUp_GeneBody.sortedByRINGBinUNT_in5kbPromoter.bed") %>%
  select(name, deepTools_group) %>%
  rename(refID = name, TimeUp = deepTools_group)

## Log2 FC for 2-8h IAA
lfc <- read.csv("../data/PRC1degron.cnRNAseq.DESeq2_results.csv") %>%
  select(refID, contains("LFC"), -contains("24h"))
names(lfc) <- gsub("LFC_", "", names(lfc))

lfc.order <- left_join(first.order, lfc, by = "refID") # the ordering based on first.order should be kept

## Expression RPKM (UNT)
rpkm <- read.csv("../data/PRC1degron.cnRNAseq.DESeq2_meanRPKM.csv") %>%
  select(refID, UNT) %>%
  rename(RNA_UNT = UNT) %>%
  mutate_if(is.numeric, log2)

rpkm.order <- left_join(first.order, rpkm, by = "refID")



# Making matrices ---------------------------------------------------------

time.mat <- first.order %>%
  column_to_rownames("refID")

lfc.mat <- lfc.order %>%
  select(-TimeUp) %>%
  column_to_rownames("refID") %>%
  as.matrix()

rpkm.mat <- rpkm.order %>%
  select(-TimeUp) %>%
  column_to_rownames("refID") %>%
  as.matrix()



# Plotting ----------------------------------------------------------------

# HM1: Basal expression (log2 RPKM in untreated cells)

## define colors for the heatmap
display.brewer.pal(9, "PiYG")
brewer.pal(9, "PiYG")
col_fun_rpkm <- colorRamp2(c(-4, 0, 4, 8), c("#C51B7D", "#F1B6DA", "#B8E186", "#4D9221"))

hm1 <- Heatmap(rpkm.mat,
  show_row_names = F,
  cluster_columns = F,
  cluster_rows = F,
  col = col_fun_rpkm,
  split = time.mat$TimeUp,
  name = "log2 (RPKM)",
  column_names_rot = 0,
  column_names_side = "top",
  column_names_centered = T,
  row_gap = unit(1, "mm"),
  height = unit(16.1, "cm"),
  width = unit(1, "cm"),
  border = T,
  use_raster = T
)
hm1


# HM2: Gene expression changes (log2 FC for treated cells)

col_fun <- colorRamp2(c(-2, 0, 2, 4), c("blue", "white", "red", "darkred"))


hm2 <- Heatmap(lfc.mat,
  show_row_names = F,
  cluster_columns = F,
  cluster_rows = F,
  col = col_fun,
  split = time.mat$TimeUp,
  name = "log2 (FC)",
  column_names_rot = 0,
  column_names_side = "top",
  column_names_centered = T,
  row_gap = unit(1, "mm"),
  height = unit(16.1, "cm"),
  width = unit(6, "cm"),
  border = T,
  use_raster = T
)
hm2


# Assembling and saving ----------------------------------------------------------

hm1 + hm2


pdf("Heatmap.DerepressedGenes.RPKM_and_LFC.pdf")
hm1 + hm2
dev.off()

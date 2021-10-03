# For Extended Data Fig. 2b, making an UpSet plot for genes with
# significantly increased expression (p-adj < 0.05, fold change > 1.5) following auxin treatment of
# PRC1deg cells at the indicated time points, or a 72 h OHT treatment of PRC1CKO cells.

library(tidyverse)
library(UpSetR)

# Extracting UP genes in degron -----------------------------------------------------

degron <- read.csv("../data/PRC1degron.cnRNAseq.DESeq2_results.csv")

up2 <- degron %>%
  filter(Reg_2h_IAA == "UP") %>%
  pull(refID)

up4 <- degron %>%
  filter(Reg_4h_IAA == "UP") %>%
  pull(refID)

up8 <- degron %>%
  filter(Reg_8h_IAA == "UP") %>%
  pull(refID)

up24 <- degron %>%
  filter(Reg_24h_IAA == "UP") %>%
  pull(refID)



# Extracting UP genes in CKO ----------------------------------------------

cko <- read.csv("../data/PRC1cko.cnRNAseq.DESeq2_results.csv")

up72 <- cko %>%
  filter(padj < 0.05) %>%
  filter(LFC_apeglm > log2(1.5)) %>%
  pull(refID)



# Making UpSet plots for all UP -------------------------------------------

all_up <- list(
  IAA_2h = up2,
  IAA_4h = up4,
  IAA_8h = up8,
  IAA_24h = up24,
  OHT_72h = up72
)

## Checking total numbers
lengths(all_up)

# sorting non-empty intersections by size, while excluding any intersection with less than 10 members for clarity.

pdf("UpSet.Degron_vs_CKO.Plus10_intersections.pdf", width = 8, height = 6, useDingbats = FALSE)
upset(fromList(all_up),
  sets = c("OHT_72h", "IAA_24h", "IAA_8h", "IAA_4h", "IAA_2h"),
  order.by = "freq",
  mb.ratio = c(0.6, 0.4),
  nintersects = 19, # manually selecting >10 members intersects
  text.scale = c(2, 1.5, 2, 1.5, 2, 1.5),
  keep.order = TRUE
)
dev.off()

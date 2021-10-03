# For Extended Data Fig. 2a, making scatterplots comparing the log2 fold changes in gene expression (cnRNA-seq)
# at different time points after auxin treatment in PRC1deg cells with
# gene expression changes in PRC1CKO (72 h OHT) cells.

library(tidyverse)
library(ggrastr)


my_theme <- theme_bw() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 14, colour = "black", margin = 4),
    plot.title = element_text(size = 12, hjust = 0.5),
    panel.border = element_rect(colour = "black", size = 0.8),
    legend.position = "none"
  )

# importing DESeq2 results
degron <- read.csv("../data/PRC1degron.cnRNAseq.DESeq2_results.csv")

cko <- read.csv("../data/PRC1cko.cnRNAseq.DESeq2_results.csv") %>%
  select(refID, padj, LFC_apeglm) %>%
  rename(
    padj_PRC1CKO = padj,
    LFC_PRC1CKO = LFC_apeglm
  )



# Correlating LFCs --------------------------------------------------------

df <- left_join(degron, cko, by = "refID")


## 2h vs CKO

corr2 <- cor.test(
  x = df$LFC_PRC1CKO,
  y = df$LFC_2h_IAA,
  method = "pearson"
)
corr2

lm2 <- summary(lm(LFC_2h_IAA ~ LFC_PRC1CKO, df))
lm2

ggplot(df, aes(LFC_PRC1CKO, LFC_2h_IAA)) +
  geom_point_rast(colour = "midnightblue", alpha = 0.3, raster.dpi = 600, size = 1) +
  stat_density2d(aes(
    fill = after_stat(level),
    alpha = after_stat(level)
  ),
  geom = "polygon",
  bins = 1000
  ) +
  scale_alpha_continuous(range = c(0, 0.8)) +
  scale_fill_gradientn(colours = colorRampPalette(c("midnightblue", "dodgerblue3", "dodgerblue2", "yellow1", "orangered", "red4"))(100)) +
  geom_smooth(method = "lm", colour = "black", linetype = "dashed", size = 0.8) +
  annotate("text",
    x = -2.5, y = 7.5, size = 6,
    label = paste0(
      "cor = ", round(corr2$estimate, digits = 2), # cor is Pearson correlation coefficient
      "\n", "R2 = ", round(lm2$r.squared, digits = 2)
    )
  ) + # R2 is the coefficient of determination for linear regression
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.6, colour = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, colour = "black") +
  scale_x_continuous(limits = c(-5, 10)) +
  scale_y_continuous(limits = c(-5, 10)) +
  xlab("log2 Fold Change 72h OHT / UNT in PRC1CKO") +
  ylab("log2 Fold Change 2h IAA / UNT in PRC1deg") +
  my_theme

ggsave("CorrelationPlot_LFC_RNA.AUX_2h_vs_PRC1cko.pdf", width = 6, height = 6)


## 4h vs CKO

corr4 <- cor.test(
  x = df$LFC_PRC1CKO,
  y = df$LFC_4h_IAA,
  method = "pearson"
)
corr4

lm4 <- summary(lm(LFC_4h_IAA ~ LFC_PRC1CKO, df))
lm4

ggplot(df, aes(LFC_PRC1CKO, LFC_4h_IAA)) +
  geom_point_rast(colour = "midnightblue", alpha = 0.3, raster.dpi = 600, size = 1) +
  stat_density2d(aes(
    fill = after_stat(level),
    alpha = after_stat(level)
  ),
  geom = "polygon",
  bins = 1000
  ) +
  scale_alpha_continuous(range = c(0, 0.8)) +
  scale_fill_gradientn(colours = colorRampPalette(c("midnightblue", "dodgerblue3", "dodgerblue2", "yellow1", "orangered", "red4"))(100)) +
  geom_smooth(method = "lm", colour = "black", linetype = "dashed", size = 0.8) +
  annotate("text",
    x = -1.5, y = 5.5, size = 6,
    label = paste0(
      "cor = ", round(corr4$estimate, digits = 2),
      "\n", "R2 = ", round(lm4$r.squared, digits = 2)
    )
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.6, colour = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, colour = "black") +
  scale_x_continuous(limits = c(-5, 10)) +
  scale_y_continuous(limits = c(-5, 10)) +
  xlab("log2 Fold Change 74h OHT / UNT in PRC1CKO") +
  ylab("log2 Fold Change 4h IAA / UNT in PRC1deg") +
  my_theme

ggsave("CorrelationPlot_LFC_RNA.AUX_4h_vs_PRC1cko.pdf", width = 6, height = 6)



## 8h vs CKO

corr8 <- cor.test(
  x = df$LFC_PRC1CKO,
  y = df$LFC_8h_IAA,
  method = "pearson"
)
corr8

lm8 <- summary(lm(LFC_8h_IAA ~ LFC_PRC1CKO, df))
lm8

ggplot(df, aes(LFC_PRC1CKO, LFC_8h_IAA)) +
  geom_point_rast(colour = "midnightblue", alpha = 0.3, raster.dpi = 600, size = 1) +
  stat_density2d(aes(
    fill = after_stat(level),
    alpha = after_stat(level)
  ),
  geom = "polygon",
  bins = 1000
  ) +
  scale_alpha_continuous(range = c(0, 0.8)) +
  scale_fill_gradientn(colours = colorRampPalette(c("midnightblue", "dodgerblue3", "dodgerblue2", "yellow1", "orangered", "red4"))(100)) +
  geom_smooth(method = "lm", colour = "black", linetype = "dashed", size = 0.8) +
  annotate("text",
    x = -1.5, y = 5.5, size = 6,
    label = paste0(
      "cor = ", round(corr8$estimate, digits = 2),
      "\n", "R2 = ", round(lm8$r.squared, digits = 2)
    )
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.6, colour = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, colour = "black") +
  scale_x_continuous(limits = c(-5, 10)) +
  scale_y_continuous(limits = c(-5, 10)) +
  xlab("log2 Fold Change 78h OHT / UNT in PRC1CKO") +
  ylab("log2 Fold Change 8h IAA / UNT in PRC1deg") +
  my_theme

ggsave("CorrelationPlot_LFC_RNA.AUX_8h_vs_PRC1cko.pdf", width = 6, height = 6)



## 24h vs CKO

corr24 <- cor.test(
  x = df$LFC_PRC1CKO,
  y = df$LFC_24h_IAA,
  method = "pearson"
)
corr24

lm24 <- summary(lm(LFC_24h_IAA ~ LFC_PRC1CKO, df))
lm24

ggplot(df, aes(LFC_PRC1CKO, LFC_24h_IAA)) +
  geom_point_rast(colour = "midnightblue", alpha = 0.3, raster.dpi = 600, size = 1) +
  stat_density2d(aes(
    fill = after_stat(level),
    alpha = after_stat(level)
  ),
  geom = "polygon",
  bins = 1000
  ) +
  scale_alpha_continuous(range = c(0, 0.8)) +
  scale_fill_gradientn(colours = colorRampPalette(c("midnightblue", "dodgerblue3", "dodgerblue2", "yellow1", "orangered", "red4"))(100)) +
  geom_smooth(method = "lm", colour = "black", linetype = "dashed", size = 0.8) +
  annotate("text",
    x = -1.5, y = 5.5, size = 6,
    label = paste0(
      "cor = ", round(corr24$estimate, digits = 2),
      "\n", "R2 = ", round(lm24$r.squared, digits = 2)
    )
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.6, colour = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, colour = "black") +
  scale_x_continuous(limits = c(-5, 10)) +
  scale_y_continuous(limits = c(-5, 10)) +
  xlab("log2 Fold Change 72h OHT / UNT in PRC1CKO") +
  ylab("log2 Fold Change 24h IAA / UNT in PRC1deg") +
  my_theme

ggsave("CorrelationPlot_LFC_RNA.AUX_24h_vs_PRC1cko.pdf", width = 6, height = 6)

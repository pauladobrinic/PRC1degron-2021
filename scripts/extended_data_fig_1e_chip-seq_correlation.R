# For Extended Data Fig. 1e, making a scatterplot showing the relationship between RING1B cChIP-seq signal
# at RING1B-bound sites in wild-type (TIR1) and PRC1deg cells before auxin treatment.

library(tidyverse)
library(ggrastr)

my_theme <- theme_bw() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
    axis.text = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 12, colour = "black", margin = 4),
    plot.title = element_text(size = 12, hjust = 0.5),
    panel.border = element_rect(colour = "black", size = 0.8),
    legend.position = "none"
  )


## import RPKM values for cChIP-seq signal
rpkm <- readRDS("../data/PRC1degron.cChIPseq_RING1B_SUZ12_MERGED.various_intervals.readCountsPerKb.RDS")

## narrow down to RING1B peaks and conditions of interest
rpkm <- rpkm %>%
  filter(intType == "RING1B peaks") %>%
  select(intPosition, RING1B_UNT, RING1B_Tir1) %>%
  rename(
    PRC1deg = RING1B_UNT,
    TIR1 = RING1B_Tir1
  )


## transform to log2 RPKM
log2rpkm <- mutate_at(rpkm, vars(PRC1deg:TIR1), list(~ log2(.)))


## calculate correlation in RING1B peaks
corr <- cor.test(
  x = log2rpkm$TIR1,
  y = log2rpkm$PRC1deg,
  method = "pearson"
)
corr


## extract linear model
lm <- summary(lm(PRC1deg ~ TIR1, log2rpkm))
lm


## Plot
ggplot(log2rpkm, aes(TIR1, PRC1deg)) +
  geom_point_rast(colour = "midnightblue", alpha = 0.3, raster.dpi = 600, size = 1) +
  stat_density2d(aes(
    fill = after_stat(level),
    alpha = after_stat(level)
  ),
  geom = "polygon",
  bins = 30
  ) +
  scale_alpha_continuous(range = c(0, 0.5)) +
  scale_fill_gradientn(colours = colorRampPalette(c("midnightblue", "dodgerblue3", "dodgerblue2", "yellow1", "orangered", "red4"))(50)) +
  geom_smooth(method = "lm", colour = "black", linetype = "dashed", size = 0.8) +
  annotate("text",
    x = 7, y = 12, size = 6,
    label = paste0(
      "cor = ", round(corr$estimate, digits = 2), # cor is Pearson correlation coefficient
      "\n", "R2 = ", round(lm$r.squared, digits = 2)  # R2 is the coefficient of determination for linear regression
    )
  ) +
  coord_cartesian(xlim = c(6, 13), ylim = c(6, 13)) +
  my_theme +
  xlab("log2 RPKM TIR1") +
  ylab("log2 RPKM PRC1deg") +
  ggtitle("RING1B")


ggsave("CorrelationPlot.RING1B_cChIP-seq_in_RING1Bpeaks.PRC1deg_vs_TIR1.pdf",
  width = 6, height = 6
)

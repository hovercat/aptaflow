#!/usr/bin/env Rscript
library("dplyr")
library("tools")
library("tidyr")
library("ggplot2")
library("showtext")
showtext_auto()

args = commandArgs(trailingOnly=TRUE)

csv_in = args[1]
out_tiff = args[2]

df_uniq_dupli <- read.table(csv_in, header=TRUE, sep="\t", quote=NULL)
df_uniq_dupli$unique <- NULL
df_uniq_dupli$duplicated <- NULL
colnames(df_uniq_dupli) <- c("Round", "Unique", "Duplicated")

# Plotting
viridis_colorpalette="plasma"

# Barplot
gg <- ggplot(
    df_uniq_dupli %>% gather(distinction, perc, -Round),
    mapping = aes(x=Round, y=perc, fill=factor(distinction))
  ) +
  geom_col() +
  labs(fill = "", colour="") +
  xlab("SELEX Round") +
  ylab("% duplication rates of individual sequences") +
  scale_fill_viridis_d(direction=-1, begin=0, end=0.8, option = viridis_colorpalette) +
  theme_bw() +
  theme(axis.text = element_text(size=11,face="bold"),
        axis.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=11, face="bold"))

ggsave(
  out_tiff,
  plot=gg,
  device=tiff(),
  width=8,
  height=8,
  units="cm"
)
dev.off()

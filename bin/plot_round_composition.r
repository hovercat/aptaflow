#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library("dplyr")
library("tools")
library("tidyr")
library("showtext")
showtext_auto()

# PLOTTING DEFAULTS
library("ggplot2")
viridis_colorpalette="plasma"
#

df_content <- read.table(args[1], sep="\t", header=TRUE)
df_content_g <- df_content %>% gather("nt", "perc", c('A','C','G','T'))

df_content_g$nt_position = df_content_g$nt_position+1
# Plot
gg_nt_content <- ggplot(data=df_content_g, mapping = aes(x=factor(nt_position, ordered=FALSE), y=perc, fill=factor(nt))) +
  geom_col() +
  labs(fill = "", colour="") +
  xlab("Nucleotide position") +
  ylab("Nucleotide proportions") +
  scale_x_discrete(breaks=c(seq(from=0, to=max(df_content_g$nt_position), by=5 ))) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_viridis_d(direction=-1, begin=0, end=0.8, option = viridis_colorpalette) +
  theme_bw() +
  theme(axis.text = element_text(size=11,face="bold"),
        axis.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"))

ggsave(
  args[2],
  plot=gg_nt_content,
  device=tiff(),
  width=8,
  height=8,
  units="cm"
)
dev.off()

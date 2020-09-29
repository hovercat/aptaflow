#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library("dplyr")
library("tools")
library("tidyr")
library("ggplot2")
library("showtext")
showtext_auto()

viridis_colorpalette="plasma"

df_content <- read.table(args[1], sep="\t", header=TRUE)
df_content_g <- df_content %>% gather("nt", "perc", c("A","C","G","T"))

# Plot
gg <- ggplot(data=df_content_g, mapping = aes(x=round, y=perc, fill=factor(nt))) +
  geom_col() +
  geom_text(
    aes(label = paste(format(perc*100, digit=2, nsmall = 1), "", sep=""), group = factor(nt)),
    position=position_stack(vjust = 0.5),
    size=4.1,
    fontface="bold",
    colour="white"
  ) +
  labs(fill = "", colour="") +
  xlab("SELEX Round") +
  ylab("Nucleotide proportions") +
  scale_fill_viridis_d(direction=-1, begin=0, end=0.8, option = viridis_colorpalette) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_bw() +
  theme(axis.text = element_text(size=11,face="bold"),
        axis.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=11, face="bold"))

ggsave(
  "selex_acgt_bargraph.tiff",
  plot=gg,
  device=tiff(),
  width=8,
  height=8,
  units="cm"
)
dev.off()



gg <- ggplot(data=df_content_g, mapping = aes(x=as.numeric(factor(round)), y=perc, fill=factor(nt))) +
  geom_area() +
  labs(fill = "", colour="") +
  xlab("SELEX Round") +
  ylab("Percentage") +
  scale_fill_viridis_d(direction=-1, begin=0, end=0.8, option = viridis_colorpalette) +
  scale_x_discrete(limits=levels(factor(df_content_g$round))) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_bw() +
  theme(axis.text = element_text(size=10,face="bold")) # pts

ggsave(
  "selex_acgt_area.tiff",
  plot=gg,
  device=tiff(),
  width=8,
  height=8,
  units="cm"
)
dev.off()

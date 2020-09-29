#!/usr/bin/env Rscript
library("dplyr")
library("tools")
library("tidyr")
library("ggplot2")
library("showtext")
showtext_auto()

args = commandArgs(trailingOnly=TRUE)

csv_in = args[1]
base = as.integer(args[2])
out_dir = args[3]


df_diversity <- read.table(csv_in, header=TRUE, sep="\t", quote=NULL)

# Plotting
viridis_colorpalette="plasma"

df_perc <- read.table(csv_in, header=TRUE, sep="\t")
df_plot <- df_perc %>%
  gather(round, perc, -exponent) %>%
  #fill(0) %>%
  mutate(exponent_b = base^exponent) %>%
  mutate(exponent_lab = paste(base^exponent, " - ", base^(exponent+1)-1, sep=""))

# Barplot
gg <- ggplot(data=df_plot, mapping = aes(x=round, y=perc, fill=factor(exponent))) +
  geom_col() +
  labs(fill="") +
  xlab("SELEX Round") +
  ylab("% total counts of individual sequences") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_viridis_d(direction=-1, begin=0, end=0.8, option = viridis_colorpalette, labels=df_plot$exponent_lab) +
  theme_bw() +
  theme(axis.text = element_text(size=11,face="bold"),
        axis.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=11, face="bold"))

ggsave(
  paste(out_dir, "round_composition_rpm_bar.tiff", sep="/"),
  plot=gg,
  device=tiff(),
  width=8,
  height=8,
  units="cm"
)
dev.off()


# Histogram
gg <- ggplot(df_plot, aes(x=as.numeric(factor(round)), y=perc, fill=factor(exponent))) +
  geom_area() +
  xlab("SELEX Round") +
  ylab("% total counts of individual sequences") +
  labs(fill = "") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_discrete(limits=levels(factor(df_plot$round))) +
  scale_fill_viridis_d(direction=-1, begin=0, end=0.8, option = viridis_colorpalette, labels=df_plot$exponent_lab) +
  theme_bw() +
  theme(axis.text = element_text(size=11,face="bold"),
        axis.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=11, face="bold"))
ggsave(
  paste(out_dir, "round_composition_rpm_hist.tiff", sep="/"),
  plot=gg,
  device=tiff(),
  width=8,
  height=8,
  units="cm"
)
dev.off()


# lineplot
gg <- ggplot(data=df_plot, mapping = aes(x=round, y=perc, group=factor(exponent), linetype=factor(exponent), shape=factor(exponent), colour=factor(exponent))) +
  geom_line(size=1) +
  geom_point() +
  labs(linetype = "Bin size", shape = "Bin size", colour = "Bin size") +
  xlab("SELEX Round") +
  ylab("% total counts of individual sequences") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_viridis_d(direction=-1, begin=0, end=0.8, option = viridis_colorpalette) +
  theme_bw() +
  theme(axis.text = element_text(size=11,face="bold"),
        axis.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=11, face="bold"))
ggsave(
  paste(out_dir, "round_composition_rpm_lines.tiff", sep="/"),
  plot=gg,
  device=tiff(),
  width=8,
  height=8,
  units="cm"
)
dev.off()

# lineplot, ylog10
gg <- ggplot(data=df_plot, mapping = aes(x=round, y=perc, group=factor(exponent), linetype=factor(exponent), shape=factor(exponent), colour=factor(exponent))) +
  geom_line(size=1) +
  geom_point() +
  labs(linetype = "Bin size", shape = "Bin size", colour = "Bin size") +
  xlab("SELEX Round") +
  ylab("% total counts of individual sequences (log10 scaled)") +
  scale_y_log10(labels = scales::percent) +
  scale_fill_viridis_d(direction=-1, begin=0, end=0.8, option = viridis_colorpalette) +
  theme_bw() +
  theme(axis.text = element_text(size=11,face="bold"),
        axis.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=11, face="bold"))
ggsave(
  paste(out_dir, "round_composition_rpm_lines_ylog10.tiff", sep="/"),
  plot=gg,
  device=tiff(),
  width=8,
  height=8,
  units="cm"
)
dev.off()

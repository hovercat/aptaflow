#!/usr/bin/env Rscript

library("argparse")
library("dplyr")
library("tools")
library("tidyr")

# Parameter validation
args_parser = ArgumentParser(
  prog="SELEXrpm",
  description=""
)
args_parser$add_argument("--in-derep-csv", "-i", help="Input file must be csv as created by SELEXderep.", required=TRUE)
args_parser$add_argument("--out-file", "-o", default="rpm.csv", help="Output file. Default: rpm.csv")
#args <- args_parser$parse_args(c("-i", "/home/boxcattu/aptaflow37/output/2020-03-02/EF05_0007/results/5_analysis/selex_derep.csv"))
args <- args_parser$parse_args()

# Data input
round_names <- c()
df_selex <-  read.csv(args$in_derep_csv, sep="\t", header=TRUE)
rownames(df_selex) <- df_selex$id
round_names <- colnames(df_selex)[c(-1:-3)]

# Set every sequence count which is 0 to NA, so the long format is more sparse.
df_selex[df_selex == 0] = NA
df_selex_long <- df_selex %>% gather("round", "count", -"id", -"seq", -"p5_seq_p3", na.rm = TRUE)
df_selex <- NULL # clear wide format from ram

df_rpm <- df_selex_long %>%
  group_by(round) %>% # calculate rpm, exponent_rpm and exponent
  mutate(
    rpm = (count/sum(count) * (10^6))
  ) %>%
  ungroup()


# write rpm to csv file
write.table(
  df_rpm,
  paste(args$out_file),
  append=FALSE,
  quote=FALSE,
  row.names = FALSE,
  sep=";",
  na="0"
)

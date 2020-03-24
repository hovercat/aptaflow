#!/usr/bin/env Rscript

library("argparse")
library("dplyr")
library("tools")
library("tidyr")


# Parameter validation
args_parser = ArgumentParser(
  prog="SELEXanalyse_duplicates",
  description="This r-script calculates the ratios between duplicate and singleton reads for every round. Input is the csv table created by SELEXderep."
)
args_parser$add_argument("--in-derep-csv", "-i", help="Input file must be csv as created by SELEXderep.", required=TRUE)
args_parser$add_argument("--out-csv", "-o", default=".", help="Name of the output csv file.", required=TRUE)
args <- args_parser$parse_args()

# Data input
round_names <- c()
df_selex <-  read.csv(args$in_derep_csv, sep="\t", header=TRUE)
rownames(df_selex) <- df_selex$id
round_names <- colnames(df_selex)[c(-1:-3)]
df_selex[df_selex == 0] = NA
df_selex <- df_selex %>% gather("round", "abundancy", -"id", -"seq", -"p5_seq_p3", na.rm = TRUE)


df_uniq_dupli <- df_selex %>%
  group_by(round) %>%
  mutate(unique = abundancy==1) %>%
  summarise(
    unique = sum(unique), 
    duplicated = n() - sum(unique),
    unique_p = sum(unique)/n(),
    duplicated_p = (n()-sum(unique))/n()
  )

write.table(
  df_uniq_dupli,
  args$out_csv,
  quote=FALSE,
  sep="\t",
  col.names=TRUE,
  row.names=FALSE
)

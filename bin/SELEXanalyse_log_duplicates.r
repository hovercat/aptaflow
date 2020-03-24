#!/usr/bin/env Rscript

library("argparse")
library("dplyr")
library("tools")
library("tidyr")

# Parameter validation
args_parser = ArgumentParser(
  prog="SELEXanalyse_log_duplicates",
  description=""
)
args_parser$add_argument("--in-rpm-csv", "-i", help="Input file must be csv as created by SELEXrpm", required=TRUE)
args_parser$add_argument("--bin-base", "-b", type="integer", default=10, help="Bins are spaced logarithmically, e.g. for bin-base of 10: [1, 10, 100, 1000, ...]. Default: 10")
args_parser$add_argument("--min-abundancy", "-m", type="integer", default=1, help="Determines how often the sequence must be present to be considered. Default: 1")
args_parser$add_argument("--out-dir", "-o", default=".", help="Output directory. Default: working directory")
#args <- args_parser$parse_args(c("-i", "rpm.csv"))
args <- args_parser$parse_args()


# Data input
df_rpm <-  read.csv(args$in_rpm_csv, sep=";", header=TRUE)
round_names <- levels(factor(df_rpm$round))
df_rpm <- df_rpm %>%
  ungroup() %>%
  group_by(round) %>% # calculate rpm, exponent_rpm and exponent
  mutate(
    rpm = (count/sum(count) * (10^6)),
    exponent_rpm = floor(log(rpm, base=args$bin_base)),
    exponent = floor(log(count, base=args$bin_base))
  ) %>%
  ungroup()

df_diversity_counted <- df_rpm %>%
  group_by(round, exponent) %>% 
  mutate(bin_size_total=sum(rpm)) %>%
  select("round", "exponent", "bin_size_total") %>% 
  distinct() %>%
  data.frame()
colnames(df_diversity_counted) = c("round", "exponent", "bin_size_total")

# Sometimes there are exponents missing, e.g. if there are multiple 0-9 seqs, no 10-99 and some 100-999 seqs.
exponents <- 0:max(df_diversity_counted$exponent)
missing_exponents <- exponents[!(exponents %in% df_diversity_counted$exponent)]
df_missing <- data.frame(
  rep(round_names, length(missing_exponents)), 
  rep(missing_exponents, each=length(round_names)), 
  rep(0, length(round_names) * length(missing_exponents))
)
colnames(df_missing) <- c("round", "exponent", "bin_size_total")
df_diversity_counted <- rbind(df_diversity_counted, df_missing, stringsAsFactors=TRUE)


write.table(
  df_diversity_counted %>% spread(round, bin_size_total, fill=0), 
  paste(args$out_dir, "rpm.log_duplicates.csv", sep="/"),
  append=FALSE,
  quote=FALSE,
  row.names = FALSE,
  sep="\t",
  na="0"
)

df_diversity_counted_percentages <- df_diversity_counted %>% mutate(bin_size_total=bin_size_total/(10^6))


write.table(
  df_diversity_counted_percentages %>% spread(round, bin_size_total, fill=0), 
  paste(args$out_dir, "rpm_percent.log_duplicates.csv", sep="/"),
  append=FALSE,
  quote=FALSE,
  row.names = FALSE,
  sep="\t",
  na="0"
)

df_diversity_unique <- df_rpm %>% 
  group_by(round, exponent) %>% 
  mutate(bin_size_unique=n()) %>%
  select("round", "exponent", "bin_size_unique") %>% 
  distinct()


write.table(
  df_diversity_unique %>% spread(round, bin_size_unique, fill=0), 
  paste(args$out_dir, "unique.log_duplicates.csv", sep="/"),
  append=FALSE,
  quote=FALSE,
  row.names = FALSE,
  sep="\t",
  na="0"
)


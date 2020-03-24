#!/usr/bin/env Rscript

library("Biostrings")
library("dplyr")
library("tools")

# # Parameter validation
# args_parser = ArgumentParser(
#   prog="SELEXntcomposition_all_rounds",
#   description="SELEXanalyse_composition counts ACGT content for every round and outputs their percentages."
# )
# args_parser.add_argument("round_compositon_csv", metavar="round_composition_csv", type=str, nargs='+', help=
#                            "CSV files in sequential order as produced by SELEXntcomposition")
# args_parser$add_argument("--out-csv", "-o", default=".", help="Output csv.", required=TRUE)
# args <- args_parser$parse_args()
# 
# 
# # Data input
# round_names <- c()
# df_selex <-  read.csv(args$in_derep_csv, sep="\t", header=TRUE)
# rownames(df_selex) <- df_selex$id
# round_names <- colnames(df_selex)[c(-1:-3)]
# df_selex[df_selex == 0] = NA
# df_selex <- df_selex %>% gather("round", "abundancy", -"id", -"seq", -"p5_seq_p3", na.rm = TRUE)


args = commandArgs(trailingOnly=TRUE)
#args[1] = 40
#args[2] = "/home/boxcattu/apta_data/resc/ef05/R07.fasta"
#args[3] = "/home/boxcattu/apta_data/resc/ef05/R08.fasta"
# args[9] = "/home/boxcattu/apta_data/resc/ef05/R09.fasta"
# args[10] = "/home/boxcattu/apta_data/resc/ef05/R10.fasta"
# args[11] = "/home/boxcattu/apta_data/resc/ef05/R11.fasta"

aptamer_n_region <- as.integer(args[1])

round_names=c()
df_selex = data.frame(round_name=character(), A=double(), C=double(), G=double(), T=double())
for (fasta_path in args[-1]) {
  round_name=basename(tools::file_path_sans_ext(fasta_path))
  round_names = c(round_names, round_name)
  cat("Processing round ")
  cat(round_name)
  cat("\n")
  
  round_letter_fq <- letterFrequency(readDNAStringSet(fasta_path), letters = c("A", "C", "G", "T"))
  df_round <- data.frame(round_letter_fq) %>%
    mutate(round = round_name) %>%
    summarise(
      round=first(round),
      A=sum(A)/(n()*aptamer_n_region),
      C=sum(C)/(n()*aptamer_n_region),
      G=sum(G)/(n()*aptamer_n_region),
      T=sum(T)/(n()*aptamer_n_region),
    )
  
  df_selex <- rbind(df_selex, df_round)
}


round_names=c()
df_selex = data.frame(round_name=character(), A=double(), C=double(), G=double(), T=double())
for (fasta_path in args[-1]) {
  round_name=basename(tools::file_path_sans_ext(fasta_path))
  round_names = c(round_names, round_name)
  cat("Processing round ")
  cat(round_name)
  cat("\n")
  
  round_letter_fq <- letterFrequency(readDNAStringSet(fasta_path), letters = c("A", "C", "G", "T"))
  df_round <- data.frame(round_letter_fq) %>%
    mutate(round = round_name) %>%
    summarise(
      round=first(round),
      A=sum(A)/(n()*aptamer_n_region),
      C=sum(C)/(n()*aptamer_n_region),
      G=sum(G)/(n()*aptamer_n_region),
      T=sum(T)/(n()*aptamer_n_region),
    )
  
  df_selex <- rbind(df_selex, df_round)
}

write.table(
  df_selex, 
  "selex.acgt.csv",
  sep="\t",
  quote=FALSE,
  col.=TRUE
)

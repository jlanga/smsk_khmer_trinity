#!/usr/bin/env Rscript

## Management of arguments
require(optparse)
option_list = list(
    make_option(
        opt_str = c("-q", "--quantification"),
        type = "character",
        default = NULL, 
        help = "tsv with the quantification from kallisto (abundancs.tsv)",
        metavar = "character"
    ),
    make_option(
        opt_str = c("-m", "--gene_to_transcript_map"),
        type = "character",
        default = NULL,
        help = "A TSV with the gene_id and transcritpt_id (the result when doing get_Trinity_gene_to_trans_map.pl < Trinity.fasta)",
        metavar = "character"
    ),
    make_option(
        opt_str = c("-e", "--expressed_ids"),
        type= "character",
        default = NULL, 
        help= "A TSV with the ids of the expressed transcripts",
        metavar="character"
    ),
    make_option(
        opt_str = c("-t", "--tpms_percent"),
        type = "numeric",
        default = 0.05,
        help= "Percentage to consider expressed at a loci (Default: 0.05 - 5%)",
        metavar= "number"
    )
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


#Main program
require(tidyverse)

# Read abundances
tpm_raw <- read_tsv(
        file= opt$quantification,
        col_names= TRUE
    ) 

# Read gene-transcript file
g2t_map <- read_tsv(
        file = opt$gene_to_transcript_map,
        col_names = c("gene_id", "target_id")
    )

# Join both files. Should be gene-transcript-tpms
tpm_full <- left_join(
    x= g2t_map,
    y= tpm_raw,
    by= "target_id"
)

# Compute the tpms per gene and the threshold
tpm_by_gene <- tpm_full %>%
    group_by(gene_id) %>%
    summarise(tpm = sum(tpm)) %>%
    mutate(low = opt$tpms_percent * tpm)

# Join again
tpm_full <- left_join(
    x= tpm_full, 
    y= tpm_by_gene
)

# Filter tpm > 0 and tpm>=low
tpm_filtered <- tpm_full %>%
    filter(tpm > 0) %>%
    filter(tpm >= low) %>%
    select(target_id)

# Write only ids
write_tsv(
    x = tpm_filtered,
    path = opt$expressed_ids,
    col_names = FALSE
)

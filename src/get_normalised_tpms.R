#!/usr/bin/env Rscript

## Management of arguments
require(optparse)
option_list = list(
    make_option(
        opt_str = c("-e", "--experimental_design"),
        type = "character",
        default = NULL, 
        help = "tsv with the experimental design. Expected column names `run_accession` and `condition` ",
        metavar = "character"
    ),
    make_option(
        opt_str = c("-s", "--sleuth_object"),
        type = "character",
        default = NULL,
        help = ".Rdata object with the prepared sleuth object (for debugging or sleuth downstream analysis)",
        metavar = "character"
    ),
    make_option(
        opt_str = c("-t", "--normalised_tpms"),
        type= "character",
        default = NULL, 
        help= "TSV with the normalised transcripts per million",
        metavar="character"
    )
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);



# Load experimental design (i.e. file - condition)
require(sleuth)
require(tidyverse)

s2c <- read_tsv(
        file = opt$experimental_design,
        col_names = TRUE
    )


# Load quant and fit model
so <- sleuth_prep(s2c, ~ tissue)

# Store results
write_tsv(so$obs_norm, path = opt$normalised_tpms)
save(so, file = opt$sleuth_object)



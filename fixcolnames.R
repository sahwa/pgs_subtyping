#!/usr/bin/env Rscript
library(data.table)
library(stringr)
setDTthreads(0L)

args = commandArgs(trailingOnly=TRUE)

setwd("/well/ckb/users/aey472/projects/pgs_subtype/data/sumstats")

x = args[1]

dat = fread(x)

fixnames = function(dat, old, new) {
	if ((old %in% colnames(dat))) colnames(dat)[which(colnames(dat) == old)] = new
	dat
}

## effect and non-effect alleles ##
dat = fixnames(dat, "EFFECT_ALLELE", "EA")
dat = fixnames(dat, "effect_allele", "EA")
dat = fixnames(dat, "Effect_A2", "EA")
dat = fixnames(dat, "Tested_Allele", "EA")
dat = fixnames(dat, "hm_effect_allele", "EA")

dat = fixnames(dat, "other_allele", "NEA")
dat = fixnames(dat, "hm_other_allele", "NEA")

## pos
dat = fixnames(dat, "position(b37)", "POS")
dat = fixnames(dat, "POS_b37", "POS")
dat = fixnames(dat, "POS_GRCh37", "POS")
dat = fixnames(dat, "Pos_b37", "POS")

## chr 
dat = fixnames(dat, "chromosome(b37)", "CHR")
dat = fixnames(dat, "#CHROM", "CHR")
## P
dat = fixnames(dat, "Fixed-effects_p-value", "P")

# SE
dat = fixnames(dat, "Fixed-effects_SE", "SE")

# BETA
dat = fixnames(dat, "Fixed-effects_beta", "BETA")

# ID
dat = fixnames(dat, "MarkerName", "MarkerID")

## try fix the ID cols 
test = dat[1:100]
rsidcols = str_detect(test, "rs[0-9]{1,}")
if (any(rsidcols)) {
  colnames(dat)[which(rsidcols)] = "rsID"
}


fwrite(dat, x, sep="\t", col.names=T, quote=F)

cat("Completed!\n")

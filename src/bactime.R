#!/usr/bin/env Rscript

library("optparse")
library(ape)
library(BactDating)
library(tidyverse)
set.seed(0)

option_list =  list(
	make_option(c("-t","--tree"), type="character",
		default=NULL, help="ML Tree", metavar="character"),
	make_option(c("-d","--date"), type="character",
		default=NULL, help="Samples Isolation date",
		metavar="character"),
	make_option(c("-l","--alnlen"), type="character",
		default=NULL, help="Alignment length",
		metavar="character"),


	make_option(c("-o","--outfile"), type="character",
		default=NULL, help="Output file for visualisation",
		metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$tree)) {
	print_help(opt_parser)
	stop("Tree file is not given.n", call.=FALSE)

}
if (is.null(opt$date)) {
	print_help(opt_parser)
	stop("Date file is not given.n", call.=FALSE)

}
if (is.null(opt$alnlen)) {
	print_help(opt_parser)
	stop("Alignment length is not given.n", call.=FALSE)

}
if (is.null(opt$outfile)) {
	print_help(opt_parser)
	stop("Outfile file is not given.n", call.=FALSE)

}

tree = read.tree(opt$tree)
tree$edge.length = tree$edge.length * as.integer(opt$alnlen)
coldate = read.csv(opt$date)
coldate = coldate[match(tree$tip.label, coldate$id)]
res = bactdate(tree,coldate, nbIts=1000000, showProgress=T)
saveRDS(file=opt$outfile)





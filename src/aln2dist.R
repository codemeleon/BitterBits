library(ape)
library(optparse)
library(stringr)
library(reshape)

option_list <- list(
  make_option(c("-a", "--alnfile"),
    type = "character", default = "~/Documents/Victor/hsdm_all.fasta",
    help = "Fasta alignment file", metavar = "character"
  ), make_option(c("-d", "--distype"),
    type = "character", default = "SNP",
    help = "Type of distance measure", metavar = "character"
  ), make_option(c("-o", "--outfile"),
    type = "character",
    default = "outputfile", help = "outputfile", metavar = "character"
  ),
  make_option(c("-s", "--sw"),
    type = "character", default = NULL,
    help = "Sliding window size", metavar = "character"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))
aln_mat <- as.matrix(read.FASTA(opt$alnfile))


if (opt$distype == "SNP") {
  dist_matrix <-
    as.matrix(dist.gene(aln_mat))
}


melted_matrix <- melt(dist_matrix)
meltted_matrix <- melted_matrix[melted_matrix$Var1 != melted_matrix$Var2]

all_dist <- TRUE

if (opt$sw == NULL) {
  write.csv(melted_matrix$value, opt$outfile, row.names = FALSE)
} else {
  if (grepl(str_sub(opt$sw, -1), "YMD", fixed = TRUE)) {
    time_val <- as.integer(str_sub(opt$sw, , -2))
    if (is.na(time_val) == TRUE) {
      stop("Time prefix are not numeric")
    } else {
      if (time_val == 0) {
        write.csv(meltted_matrix$value, opt$outfile, row.names = FALSE)
      } else {

      }
      # TODO: Merge time
      # TODO: Convert time in days, month and year
    }
  } else {
    stop("Given time scale is not Year, Month or Day")
  }
}





if (all_dist == TRUE) {
  write.csv(meltted_matrix$value, opt$outfile, row.names = FALSE)
  # Write about output file
}

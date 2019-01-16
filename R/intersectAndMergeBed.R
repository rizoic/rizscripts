#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(rizlib))
suppressPackageStartupMessages(library(futile.logger))
suppressPackageStartupMessages(library(GenomicRanges))

option_list = list(
  make_option(c("-b", "--bedfiles"), action="store", default=NA, type='character',
              help="A comma seperated list of bed files to intersect and merge"),
  make_option(c("-m", "--minOverlap"), action="store", default=NA, type='character',
              help="Number of files the peak should be in. Or proportion of files it should be in [default all]"),
  make_option(c("-o", "--outfile"), action="store", default=NA, type='character',
              help="Output bed file with intersect and merge results [default STDOUT]")  
)

opt = parse_args(OptionParser(option_list=option_list))

# Split the comman seperated bed file list
bed.files <- unlist(strsplit(x = opt$bedfiles, split = ",", fixed = TRUE))

# If minoverlap specfied use it else set to all the bed files
if(is.na(opt$minOverlap)){
  minOverlap = length(bed.files)
} else{
  minOverlap = opt$minOverlap
}

# Get the intersecting and merged peaks
out.bed.grange <- rizlib::intersectAndExtendBed(bedFiles = bed.files,
                                minOverlap = as.numeric(minOverlap))

# Convert to a df for easier handling
out.bed.df <- data.frame(chrom = GenomicRanges::seqnames(out.bed.grange),
                         start = GenomicRanges::start(out.bed.grange),
                         end = GenomicRanges::end(out.bed.grange))

# If no output file specified print to STDOUT
if(is.na(opt$outfile)){
  write.table(x=out.bed.df, file="", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
}else{
  write.table(x=out.bed.df, file=opt$outfile, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
}
#!/usr/bin/env Rscript
library(optparse)
library(ATACseqQC)
option_list = list(
    make_option(c("-b", "--bam_file_path"), type="character", default=NULL, 
                help="input bam file, must be indexed first and duplicates marked", metavar="character"),
    make_option(c("-o", "--output"), type="character", default = NULL,
                help="output directory for plots", metavar="character"),
    make_option(c("-p", "--prefix"), type="character", default = NULL,
                help="prefix for plotting files", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

bam_file_name = bam_file_name = gsub(".bam", "", basename(opt$bam_file_path))

#message("Calculating read duplcation frequency")
rdf = readsDupFreq(opt$bam_file_path)

pdf(file.path(opt$output, paste0(opt$prefix, '_QC_plots.pdf')))
try(hist(as.vector(rep(rdf[, 'Var1'] , rdf[, 'Freq'])), 
    breaks = 200, 
    ylim = c(0, 300), 
    col = 'darkolivegreen3', 
    xlab = 'Number of Duplicates', 
    ylab = 'Frequency', 
    main = 'Reads duplcation frequency'), silent = TRUE)       #plot 1
try(estimateLibComplexity(rdf), silent = TRUE)                     #plot 2
try(fragSizeDist(opt$bam_file_path, bam_file_name), silent = TRUE)     #plot 3
dev.off()


#!/usr/bin/env Rscript
library(consensusSeekeR)
library(argparse)
library(rtracklayer)
parser <- ArgumentParser()
parser$add_argument('-n', '--narrowPeaks', nargs='+',help='list of narrowpeak files')
parser$add_argument('-g', '--genomeInfo',help='path to genome info tsv (can be generated using samtools idxstats)')
parser$add_argument('-o', '--out',help='path to write output')
parser$add_argument('-p', '--ncores', help='number of cores to use')
args <- parser$parse_args()
genome_info = read.table(file = args$genomeInfo col.names = c('seq_name', 'seq_len', 'mapped_reads', 'unmapped_reads'))
genome_info = genome_info[-nrow(genome_info),]
genome_info = Seqinfo(seqnames = genome_info[['seq_name']], 
                        seqlengths = genome_info[['seq_len']], 
                        isCircular = rep(FALSE, times = nrow(genome_info)))
peaks = c()
summits = c()
for(file_narrowPeak in args$narrowPeaks) {
    result = readNarrowPeakFile(file_narrowPeak, extractRegions = TRUE, extractPeaks = TRUE)
    peaks[basename(file_narrowPeak)] = result$narrowPeak
    summits[basename(file_narrowPeak)] = result$peak
}
consensus_regions = findConsensusPeakRegions(narrowPeaks, 
                                                peaks, 
                                                chrInfo, 
                                                extendingSize = 75,
                                                expandToFitPeakRegion = FALSE, 
                                                shrinkToFitPeakRegion = TRUE,
                                                minNbrExp = 1L, 
                                                nbrThreads = parser$ncorse)
export.bed(consensus_regions$consensusRanges,file=args$out)
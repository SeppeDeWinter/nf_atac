#!/usr/bin/env Rscript
library(consensusSeekeR)
library(argparse)
library(rtracklayer)
parser <- ArgumentParser()
parser$add_argument('-n', '--narrowPeaks', nargs='+',help='list of narrowpeak files')
parser$add_argument('-g', '--genomeInfo',help='path to genome info tsv (can be generated using samtools idxstats)')
parser$add_argument('-o', '--out',help='path to write output')
parser$add_argument('-p', '--ncores', help='number of cores to use', type = "integer")
args <- parser$parse_args()
print(args)
genome_info = read.table(file = args$genomeInfo, col.names = c('seq_name', 'seq_len', 'mapped_reads', 'unmapped_reads'))
genome_info = genome_info[-nrow(genome_info),]
genome_info = Seqinfo(seqnames = genome_info[['seq_name']], 
                        seqlengths = genome_info[['seq_len']], 
                        isCircular = rep(FALSE, times = nrow(genome_info)))
narrowPeaks = GRanges()
peaks = GRanges()
for(file_narrowPeak in args$narrowPeaks) {
    result = readNarrowPeakFile(file_narrowPeak, extractRegions = TRUE, extractPeaks = TRUE)
    
    narrowPeak = result$narrowPeak
    names(narrowPeak) = rep(basename(file_narrowPeak), length(narrowPeak))
    narrowPeaks = suppressWarnings(c(narrowPeaks, narrowPeak))

    peak = result$peak
    names(peak) = rep(basename(file_narrowPeak), length(peak))
    peaks = suppressWarnings(c(peaks, peak))
}

consensus_regions = findConsensusPeakRegions(narrowPeaks, 
                                                peaks, 
                                                genome_info, 
                                                extendingSize = 75,
                                                expandToFitPeakRegion = TRUE, 
                                                shrinkToFitPeakRegion = TRUE,
                                                minNbrExp = 1L, 
                                                nbrThreads = args$ncores)
gr = consensus_regions$consensusRanges
df <- data.frame(seqnames=seqnames(gr),
  starts=start(gr)-1,
  ends=end(gr),
  names=c(rep(".", length(gr))),
  scores=c(rep(0, length(gr))),
  strands=strand(gr))

write.table(df, file=args$out, quote=F, sep="\t", row.names=F, col.names=F)

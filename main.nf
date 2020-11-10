/**
==modules===
conda env nf_ATAC
Subread
**/

bams = Channel.fromPath(params.raw_bam_files_path)

process fixmate_and_markdup {
    label 'default'
    input:
        path RAW_BAM_FILE from bams
    output:
        path "MARKDUP-${RAW_BAM_FILE.baseName}.bam" into markdupRawBam_ch1
        path "MARKDUP-${RAW_BAM_FILE.baseName}.bam" into markdupRawBam_ch2
        path "MARKDUP-${RAW_BAM_FILE.baseName}.bam" into markdupRawBam_ch3
        path "MARKDUP-${RAW_BAM_FILE.baseName}.bam" into markdupRawBam_ch4
    """
    samtools collate -o collate_${RAW_BAM_FILE.baseName}.bam ${RAW_BAM_FILE} -@ ${task.cpus} tmp_${RAW_BAM_FILE.baseName}.bam
    samtools fixmate -m collate_${RAW_BAM_FILE.baseName}.bam -@ ${task.cpus} - | samtools sort -O bam -@ ${task.cpus} - | samtools markdup -@ ${task.cpus} - MARKDUP-${RAW_BAM_FILE.baseName}.bam
    """
} 

process pre_filter_QC {
    label 'default'
    input:
        path MARKDUP_RAW_BAM_FILE from markdupRawBam_ch1
        path QC_DIR from params.qc_dir
        path src from params.src_folder
    shell:
    '''
    samtools index -@ !{task.cpus} !{MARKDUP_RAW_BAM_FILE}
    sample_name=`echo !{MARKDUP_RAW_BAM_FILE.baseName} | awk -F'-' '{print $NF}'`
    Rscript !{src}/ATAC_QC.R \
            -b !{MARKDUP_RAW_BAM_FILE} \
            -o !{QC_DIR} \
            -p pre_filter_${sample_name}
    '''
}

process pre_filter_fingerplot {
    label 'default'
    input:
        path BAM_FILE from markdupRawBam_ch2.collect()
        path QC_DIR from params.qc_dir
    shell:
    '''
    for b in !{BAM_FILE}
    do
        samtools index -@ !{task.cpus} $b
    done
    plotFingerprint -b !{BAM_FILE} -p !{task.cpus} -plot !{QC_DIR}/pre_filter_FingerPrintPlot.pdf
    '''
    
}

process filter {
    label 'default'
    input:
        path MARKDUP_RAW_BAM_FILE from markdupRawBam_ch3
        path OUT_DIR from params.out_dir
        val MAPQ_THRESH from params.mapQ_thr
        val MITO_NAME from params.mito_name
    output:
        path "${OUT_DIR}/noMito_Q${MAPQ_THRESH}_FILTERED-${MARKDUP_RAW_BAM_FILE.baseName}.bam" into filteredBam_ch1
        path "${OUT_DIR}/noMito_Q${MAPQ_THRESH}_FILTERED-${MARKDUP_RAW_BAM_FILE.baseName}.bam" into filteredBam_ch2
        path "${OUT_DIR}/noMito_Q${MAPQ_THRESH}_FILTERED-${MARKDUP_RAW_BAM_FILE.baseName}.bam" into filteredBam_ch3
        path "${OUT_DIR}/noMito_Q${MAPQ_THRESH}_FILTERED-${MARKDUP_RAW_BAM_FILE.baseName}.bam" into filteredBam_ch4
    shell:
    '''
    samtools index !{MARKDUP_RAW_BAM_FILE} 
    samtools idxstats !{MARKDUP_RAW_BAM_FILE} | cut -f 1 | grep -v !{MITO_NAME} | xargs samtools view -q !{MAPQ_THRESH} -u !{MARKDUP_RAW_BAM_FILE} | samtools sort /dev/stdin -@ !{task.cpus} -o !{OUT_DIR}/noMito_Q!{MAPQ_THRESH}_FILTERED-!{MARKDUP_RAW_BAM_FILE.baseName}.bam
    '''

}

process call_peaks_filtered {
    label 'default'
    input:
        path bam from filteredBam_ch3
        path OUT_DIR from params.out_dir
        val GENOME_SIZE from params.macs_effective_genome_size
    output:
        path "${OUT_DIR}/filtered_${bam.baseName.tokenize('.').get(0)}_peaks.narrowPeak" into bed_ch
    shell:
    '''
    macs2 callpeak -t !{bam} -g !{GENOME_SIZE} --outdir !{OUT_DIR} -n filtered_!{bam.baseName.tokenize('.').get(0)} --nomodel
    '''
}

process make_consensus_peak_set {
    label 'default'
    input:
        path beds from bed_ch.collect()
        path GENOME_INFO from params.genome_chr_sizes
        path src from params.src_folder
        path OUT_DIR from params.out_dir
    output:
        path "${OUT_DIR}/consensus.bed" into consensus_bed_ch
    shell:
    '''
    Rscript !{src}/make_consensus_peaks.R \
            -n !{beds} \
            -g !{GENOME_INFO} \
            -o !{OUT_DIR}/consensus.bed \
            -p !{task.cpus}
    '''
}

process count_reads_in_consensus_peaks {
    label 'default'
    input:
        path consensus_bed from consensus_bed_ch
        path src from params.src_folder
        path bam from filteredBam_ch4
    output:
        tuple path(bam), path("${bam.baseName.tokenize('.').get(0)}_RIP.summary") into bam_rip_ch
    shell:
    '''
    GFF_FILE=!{consensus_bed.baseName.tokenize('.').get(0)}.gff
    !{src}/bedToGFF.py !{consensus_bed} ${GFF_FILE}
    featureCounts !{bam} \
                   -a ${GFF_FILE} \
                   -T !{task.cpus} \
                   -f -O \
                   -g region \
                   -t Region \
                   -o !{bam.baseName.tokenize('.').get(0)}_RIP
    '''
}

process make_FRIP_norm_BIGWIGS {
    label 'default'
    input:
        tuple path(bam), path(count_summary) from bam_rip_ch
        path OUT_DIR from params.out_dir
        val GENOME_SIZE from params.macs_effective_genome_size
    output:
        path "${OUT_DIR}/${bam.baseName.tokenize('.').get(0)}.bw" into bigwig_ch
    shell:
    '''
    RIP=`cat !{count_summary} | awk 'NR==2' | cut -f 2`
    SCALING_FACTOR=$(echo "scale=9; 1/($RIP/1000000)" | bc)
    for bam in !{bam}
    do
        samtools index $bam
    done
    bamCoverage \
        --bam !{bam} \
        --scaleFactor ${SCALING_FACTOR} \
        --normalizeUsing None \
        --binSize 1 \
        -p !{task.cpus} \
        --effectiveGenomeSize !{GENOME_SIZE} \
        -o !{OUT_DIR}/!{bam.baseName.tokenize('.').get(0)}.bw
    '''
}

process TSS_plot {
    label 'default'
    input:
        path QC_DIR from params.qc_dir
        path TSS_bed from params.TSS_bed
        path bigwigs from bigwig_ch.collect()
    shell:
    '''
    computeMatrix reference-point \
                  -S !{bigwigs} \
                  -R !{TSS_bed} \
                  -o bw_TSS_matrix.gz \
                  --referencePoint TSS \
                  -b 3000 -a 3000 \
                  -bs 5 \
                  --smartLabels \
                  -p !{task.cpus}
    plotProfile -m bw_TSS_matrix.gz \
                -o !{QC_DIR}/TSS_plot.pdf 
    '''
}
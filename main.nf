/**
==modules===
samtools
conda env nf_ATAC
deepTools
macs2
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
        path "${OUT_DIR}/noMito_Q${MAPQ_THRESH}_FILTERED-${MARKDUP_RAW_BAM_FILE.baseName.tokenize('-').get(1)}.bam" into filteredBam_ch1
        path "${OUT_DIR}/noMito_Q${MAPQ_THRESH}_FILTERED-${MARKDUP_RAW_BAM_FILE.baseName.tokenize('-').get(1)}.bam" into filteredBam_ch2
        path "${OUT_DIR}/noMito_Q${MAPQ_THRESH}_FILTERED-${MARKDUP_RAW_BAM_FILE.baseName.tokenize('-').get(1)}.bam" into filteredBam_ch3
    shell:
    '''
    samtools index !{MARKDUP_RAW_BAM_FILE} 
    samtools idxstats !{MARKDUP_RAW_BAM_FILE} | cut -f 1 | grep -v !{MITO_NAME} | xargs samtools view -q !{MAPQ_THRESH} -u !{MARKDUP_RAW_BAM_FILE} | samtools sort /dev/stdin -@ !{task.cpus} -o !{OUT_DIR}/noMito_Q!{MAPQ_THRESH}_FILTERED-!{MARKDUP_RAW_BAM_FILE.baseName.tokenize('-').get(1)}.bam
    '''

}

process call_peaks_filtered {
    label 'default'
    input:
        tuple val(REP), file(BAMS) from filteredBam_ch3.map { p -> [p.baseName.tokenize('-').get(1).tokenize('_').get(0), p]}.groupTuple()
        path OUT_DIR from params.out_dir
        val Q_VAL from params.macs_q
        val SHIFT from params.macs_shift
        val GENOME_SIZE from params.macs_effective_genome_size
        val EXTENSION_SIZE from params.macs_extension_size
    output:
        path "filtered_${REP}_q${Q_VAL}.narrowPeak" into filteredBed_ch
    shell:
    '''
    macs2 callpeak -t !{BAMS} \
    -q !{Q_VAL} \
    -n filtered_!{REP}_q!{Q_VAL} \
    --outdir !{OUT_DIR} \
    -f BAMPE \
    -g !{GENOME_SIZE} \
    --nomodel \
    --shift !{SHIFT} \
    --extsize !{EXTENSION_SIZE} \
    --keep-dup all \
    --call-summits
    '''
}


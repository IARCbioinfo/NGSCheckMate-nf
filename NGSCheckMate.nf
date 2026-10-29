#!/usr/bin/env nextflow

// Copyright (C) 2026 IARC/WHO
// This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details <http://www.gnu.org/licenses/>.

nextflow.enable.dsl=2

// ---------------------------
// PARAMETERS
// ---------------------------

params.help          = null
params.input_folder  = null
params.input         = null
params.input_file    = null
params.output_folder = "."
params.ref           = null
params.bed           = null
params.index_ext     = ".cram.crai"
params.ext           = ".cram"
params.bam           = null
params.mem           = 16
params.cpu           = 4

// Header for the IARC tools - logo generated using the following page : http://patorjk.com/software/taag  (ANSI logo generator)
def IARC_Header() {
    return """
#################################################################################
# ██╗ █████╗ ██████╗  ██████╗██████╗ ██╗ ██████╗ ██╗███╗   ██╗███████╗ ██████╗  #
# ██║██╔══██╗██╔══██╗██╔════╝██╔══██╗██║██╔═══██╗██║████╗  ██║██╔════╝██╔═══██╗ #
# ██║███████║██████╔╝██║     ██████╔╝██║██║   ██║██║██╔██╗ ██║█████╗  ██║   ██║ #
# ██║██╔══██║██╔══██╗██║     ██╔══██╗██║██║   ██║██║██║╚██╗██║██╔══╝  ██║   ██║ #
# ██║██║  ██║██║  ██║╚██████╗██████╔╝██║╚██████╔╝██║██║ ╚████║██║     ╚██████╔╝ #
# ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═════╝ ╚═╝ ╚═════╝ ╚═╝╚═╝  ╚═══╝╚═╝      ╚═════╝  #
# Nextflow pipelines for cancer genomics.########################################
"""
}

// ---------------------------
// PROCESSES
// ---------------------------

process BCFTOOLS_CALLING {

    tag "${ID}${suffix}"

    cpus   params.cpu
    memory "${params.mem}G"

    publishDir "${params.output_folder}/vcfs/", mode: 'copy'

    input:
        tuple val(ID), val(suffix), path(bam), path(bai)
        path genome
        path genome_fai
        path bed

    output:
        path("*.vcf"),            emit: vcf
        path("*.vcf.gz*"),        emit: vcfgz
        path("input_plots*.tsv"), emit: plot_inputs

    shell:
    cpus_mpileup = params.cpu.intdiv(2)
    cpus_call    = params.cpu.intdiv(2)
    file_tag     = bam.baseName
    '''
    bcftools mpileup --threads !{cpus_mpileup} --max-depth 5000 -Ou -I -R !{bed} -f !{genome} !{bam} | \
        bcftools call --threads !{cpus_call} -c -o !{file_tag}_allSM.vcf
    for sample in $(bcftools query -l !{file_tag}_allSM.vcf); do
        bcftools view -Ou -s ${sample} !{file_tag}_allSM.vcf | bcftools sort -Ou | bcftools norm -d none -O v -o ${sample}!{suffix}.vcf
        bcftools view -Oz -o ${sample}!{suffix}.vcf.gz ${sample}!{suffix}.vcf
        tabix -p vcf ${sample}!{suffix}.vcf.gz
        echo "${sample}!{suffix}.vcf\t!{ID}\t${sample}!{suffix}" >> input_plots_${sample}!{suffix}.tsv
    done
    rm !{file_tag}_allSM.vcf
    '''
}

process NCM_RUN {

    cpus   1
    memory "${params.mem}G"

    publishDir params.output_folder, mode: 'copy'

    input:
        path(vcf)
        path(genome)
        path(bed)

    output:
        path("NCM_output"), emit: ncm
        path("listVCF"),    emit: vcflist

    script:
    """
    cp -r \$NCM_HOME .
    export NCM_HOME="\$PWD/NGSCheckMate"
    echo "REF=${params.ref}" > \$NCM_HOME/ncm.conf
    echo "SAMTOOLS=samtools" >> \$NCM_HOME/ncm.conf
    echo "BCFTOOLS=bcftools" >> \$NCM_HOME/ncm.conf
    ls \$PWD/*.vcf > listVCF

    mkdir NCM_output

    python \$NCM_HOME/ncm.py -V -l listVCF -bed ${bed} -O ./NCM_output
    """
}

process NCM_GRAPHS {

    cpus   1
    memory '2 GB'

    publishDir "${params.output_folder}/NCM_output", mode: 'copy'

    input:
        path("NCM_output")
        path(infile)
        path(graphfiles)

    output:
        path("*.xgmml"), emit: graphs

    shell:
    '''
    cat !{infile} > input_plots.tsv
    Rscript !{projectDir}/bin/plots.R input_plots.tsv
    '''
}

// ---------------------------
// WORKFLOW
// ---------------------------

workflow {

    log.info IARC_Header()
    log.info ""
    log.info "-------------------------------------------------------------------------"
    log.info "    NGSCheckMate-nf : Check matching of samples         "
    log.info "-------------------------------------------------------------------------"
    log.info "Copyright (C) IARC/WHO"
    log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
    log.info "This is free software, and you are welcome to redistribute it under certain conditions; see LICENSE for details."
    log.info "-------------------------------------------------------------------------"
    log.info ""

    // --------------------------------------------------
    // HELP
    // --------------------------------------------------

    if (params.help) {
        log.info "---------------------------------------------------------------------"
        log.info "  USAGE                                                 "
        log.info "---------------------------------------------------------------------"
        log.info ""
        log.info "nextflow run iarcbioinfo/NGSCheckMate-nf [OPTIONS]"
        log.info ""
        log.info "Mandatory arguments:"
        log.info "--input_folder      FOLDER        Folder with BAM/CRAM files"
        log.info "--ref               FASTA FILE    Reference FASTA file"
        log.info ""
        log.info "Optional arguments:"
        log.info "--input             FILES         List of BAM/CRAM files (between quotes)"
        log.info '--input_file        FILE          Tab-separated file with 3 columns:'
        log.info '                                  ID, suffix, BAM/CRAM path.'
        log.info "--bed               BED FILE      Selected SNPs file (default: SNP_GRCh38.bed)"
        log.info "--output_folder     FOLDER        Output folder for NCM results [default: .]"
        log.info "--bam               FLAG          Use .bam / .bam.bai instead of .cram / .cram.crai"
        log.info "--ext               STRING        Custom alignment file extension [default: .cram]"
        log.info "--index_ext         STRING        Custom index file extension [default: .cram.crai]"
        log.info "--mem               INTEGER       Memory in GB [default: 16]"
        log.info "--cpu               INTEGER       Number of threads [default: 4]"
        exit 0
    }

    // --------------------------------------------------
    // PARAMETER CHECKS
    // --------------------------------------------------

    if (!params.ref) {
        exit 1, "Reference FASTA file not provided (--ref)."
    }

    if (!params.input_folder && !params.input_file && !params.input) {
        exit 1, "No input provided: use --input_folder, --input_file, or --input."
    }

    // --------------------------------------------------
    // LOG PARAMETERS
    // --------------------------------------------------

    log.info "input_folder  = ${params.input_folder}"
    log.info "input_file    = ${params.input_file}"
    log.info "input         = ${params.input}"
    log.info "ref           = ${params.ref}"
    log.info "output_folder = ${params.output_folder}"
    log.info "bed           = ${params.bed}"
    log.info "bam           = ${params.bam}"
    log.info "mem           = ${params.mem}"
    log.info "cpu           = ${params.cpu}"
    log.info ""

    // --------------------------------------------------
    // RESOLVE EXTENSIONS  (--bam flag overrides defaults)
    // --------------------------------------------------

    def ext       = params.bam ? ".bam"     : params.ext
    def index_ext = params.bam ? ".bam.bai" : params.index_ext

    // --------------------------------------------------
    // BUILD INPUT CHANNEL
    // --------------------------------------------------

    if (params.input_folder) {
        log.info "Input mode: folder"
        inputbams = Channel
            .fromFilePairs("${params.input_folder}/*{${ext},${index_ext}}")
            .map { id, files -> tuple(id, "", files[0], files[1]) }

    } else if (params.input_file) {
        log.info "Input mode: TSV file"
        inputbams = Channel
            .fromPath(params.input_file)
            .splitCsv(header: true, sep: '\t', strip: true)
            .map { row -> tuple(row.ID, row.suffix, file(row.bam), file(row.bam + index_ext)) }

    } else {
        log.info "Input mode: file list"
        inputbams = Channel
            .fromPath(params.input)
            .map { f -> tuple(f.baseName, "", f, f.parent / f.baseName + index_ext) }
    }

    // --------------------------------------------------
    // REFERENCE AND BED CHANNELS
    // --------------------------------------------------

    ch_ref     = Channel.value(file(params.ref))
    ch_ref_fai = Channel.value(file(params.ref + '.fai'))
    ch_bed     = params.bed
                     ? Channel.value(file(params.bed))
                     : Channel.value(file("${projectDir}/SNP_GRCh38.bed"))

    ch_graphfiles = Channel.fromPath("${projectDir}/bin/graph/*").collect()

    // --------------------------------------------------
    // RUN PROCESSES
    // --------------------------------------------------

    bcftools_out = BCFTOOLS_CALLING(inputbams, ch_ref, ch_ref_fai, ch_bed)

    // Print collected VCFs (replaces .subscribe{} from DSL1)
    bcftools_out.vcf.collect().view { files -> "VCFs produced: ${files}" }

    ncm_out = NCM_RUN(
        bcftools_out.vcf.collect(),
        ch_ref,
        ch_bed
    )

    NCM_GRAPHS(
        ncm_out.ncm,
        bcftools_out.plot_inputs.collect(),
        ch_graphfiles
    )
}

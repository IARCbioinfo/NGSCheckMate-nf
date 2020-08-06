#! /usr/bin/env nextflow

// Copyright (C) 2018 IARC/WHO

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

params.help = null

// Parameters Init
params.input_folder  = null
params.input         = null
params.input_file    = null
params.output_folder = "."
params.ref           = null
params.bed           = null
params.bai_ext       = ".bam.bai"
params.mem           = 16
params.cpu           = 4

log.info ""
log.info "-------------------------------------------------------------------------"
log.info "    NGSCheckMate-nf v1.1: Check matching of samples         "
log.info "-------------------------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "-------------------------------------------------------------------------"
log.info ""

if (params.help)
{
    log.info "---------------------------------------------------------------------"
    log.info "  USAGE                                                 "
    log.info "---------------------------------------------------------------------"
    log.info ""
    log.info "nextflow run iarcbioinfo/NGSCheckMate-nf [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--input_folder           FOLDER             Folder with BAM files"
    log.info "--ref                    FASTA FILE         Reference FASTA file"
    log.info ""
    log.info "Optional arguments:"
    log.info "--input                  BAM FILES          List of BAM files (between quotes)"
    log.info '--input_file             STRING             Input file (comma-separated) with 3 columns:'
    log.info '                                            ID (individual ID), suffix (suffix for sample names; e.g. RNA),'
    log.info '                                            and bam (path to bam file).'
    log.info "--bed                    BED FILE           Selected SNPs file (default: SNP_GRCh38.bed from the workflow's directory)"
    log.info "--output_folder          FOLDER             Output for NCM results"
    log.info "--bai_ext                STRING             Extenstion of bai files (default: .bam.bai)"
    log.info "--mem                    INTEGER            Memory (in GB)"
    log.info "--cpu                    INTEGER            Number of threads for germline calling"
    exit 0
}else{
  log.info "input_folder  = ${params.input_folder}"
  log.info "input_file    = ${params.input_file}"
  log.info "input         = ${params.input}"
  log.info "ref           = ${params.ref}"
  log.info "output_folder = ${params.output_folder}"
  log.info "bed           = ${params.bed}"
  log.info "bai_ext       = ${params.bai_ext}"
  log.info "mem           = ${params.mem}"
  log.info "cpu           = ${params.cpu}"
  log.info "help          = ${params.help}"
}


//
// Parse Input Parameters
//
if(params.input_folder){
	println "folder input"
	bam_ch = Channel.fromFilePairs("${params.input_folder}/*{.bam,$params.bai_ext}")
                         .map { row -> tuple(row[0],"",row[1][0], row[1][1]) }
}else{
    if(params.input_file){
        println "TSV file list input"
        bam_ch = Channel.fromPath("${params.input_file}")
			            .splitCsv(header: true, sep: '\t', strip: true)
			            .map { row -> [row.ID , row.suffix , file(row.bam), file(row.bam+'.bai') ] }
    }else{
	    println "file input"
	    if(params.input){
		    bam_ch = Channel.fromPath(params.input)
			                .map { input -> tuple(input.baseName, "", input, input.parent / input.baseName + '.bai') }
	    }
    }
}
ref       = file(params.ref)
if(params.bed){
    bed   = file(params.bed)
}else{
    bed   = file("${baseDir}/SNP_GRCh38.bed")
}
ncm_graphfiles = Channel.fromPath("$baseDir/bin/graph/*")

//
// Process Calling on SNP regions
//
process BCFTOOLS_calling{
    tag "$file_tag"

    cpus params.cpu
    memory params.mem+'G'

    input:
    file genome from ref 
    set ID, suffix, file(bam), file(bai) from bam_ch
    file bed

    output:
    file("*.vcf") into vcf_ch
    file("*.vcf.gz*") into vcfgz_ch
    file("input_plots*.tsv") into plot_inputs

    publishDir params.output_folder+"/vcfs/", mode: 'copy'

    shell:
    cpus_mpileup = params.cpu.intdiv(2)
    cpus_call = params.cpu.intdiv(2)
    file_tag = bam.name.replace(".bam","")
	'''
    samtools faidx !{genome}
    bcftools mpileup --threads !{cpus_mpileup} --max-depth 5000 -Ou -I -R !{bed} -f !{genome} !{bam} | bcftools call --threads !{cpus_call} -c -o !{file_tag}_allSM.vcf
    for sample in `bcftools query -l !{file_tag}_allSM.vcf`; do
        bcftools view -Ou -s ${sample} !{file_tag}_allSM.vcf | bcftools sort -Ou | bcftools norm -d none -O v -o ${sample}!{suffix}.vcf
        bcftools view -Oz -o ${sample}!{suffix}.vcf.gz ${sample}!{suffix}.vcf
        tabix -p vcf ${sample}!{suffix}.vcf.gz
        echo "${sample}!{suffix}.vcf\t!{ID}\t${sample}!{suffix}" >> input_plots_${sample}!{suffix}.tsv
    done
    rm !{file_tag}_allSM.vcf
    '''
}

vcf_ch4print = Channel.empty()
vcf_ch2 = Channel.empty()

vcf_ch.into{ vcf_ch2; vcf_ch4print}

vcf_ch4print.collect()
            .subscribe{ row -> println "${row}" }


process NCM_run {
    cpus 1
    memory params.mem+'G'

    publishDir params.output_folder, mode: 'copy'
    
    input:
    file (vcf) from vcf_ch2.collect()
    file genome from ref 

	output:
    file ("NCM_output") into ncm_ch
    file ("listVCF") into vcflist
	
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

process NCM_graphs {
    cpus 1
    memory '2 GB'

    publishDir "${params.output_folder}/NCM_output", mode: 'copy'
    
    input:
    file ("NCM_output") from ncm_ch
    file infile from plot_inputs.collect()
    file graphfiles from ncm_graphfiles.collect()

	output:
    file ("*.xgmml") into graphs

    shell:
	'''
    cat !{infile} > input_plots.tsv
    Rscript !{baseDir}/bin/plots.R input_plots.tsv
	'''
}
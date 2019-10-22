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
params.output_folder = "."
params.ref           = null
params.bed           = null
params.bai_ext       = ".bam.bai"
params.NCM_labelfile = 'NO_FILE'

log.info ""
log.info "-------------------------------------------------------------------------"
log.info "    NGSCheckMate-nf v1: Test cohort for duplicated samples       "
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
    log.info "--input                  BAM FILES             List of BAM files (between quotes)"
    log.info "--output_folder          FOLDER         Output for NCM results"
    log.info "--ref                    FASTA FILE            Reference FASTA file"
    log.info "--bed                    BED FILE              Selected SNPs file"
    exit 0
}else{
  log.info "input_folder=${params.input_folder}"
  log.info "input=${params.input}"
  log.info "ref=${params.ref}"
  log.info "output_folder=${params.output_folder}"
  log.info "bed=${params.bed}"
  log.info "help=${params.help}"
}


//
// Parse Input Parameters
//
if(params.input_folder){
	println "folder input"
	bam_ch = Channel.fromFilePairs("${params.input_folder}/*{.bam,$params.bai_ext}")
                         .map { row -> tuple(row[0],row[1][0], row[1][1]) }

	//bam_ch4print = Channel.fromFilePairs("${params.input_folder}/*{.bam,$params.bai_ext}")
    //                      .map { row -> tuple(row[0],row[1][0], row[1][1]) }
	//		  .subscribe { row -> println "${row}" }
}else{
	println "file input"
	if(params.input){
		bam_ch = Channel.fromPath(params.input)
			.map { input -> tuple(input.baseName, input, input.parent / input.baseName + '.bai') }
	}
}
output    = file(params.output_folder)
ref       = file(params.ref)
bed       = file(params.bed)
labelfile = file(params.NCM_labelfile)

//
// Process Calling on SNP regions
//
process BCFTOOLS_calling{
    tag "$sampleID"

    cpus 1
    memory '16 GB'
    time { (2.hour + (2.hour * task.attempt)) }

    errorStrategy 'retry'
    maxRetries 3

    input:
    file genome from ref 
    set sampleID, file(bam), file(bai) from bam_ch
    file bed

    output:
        set sampleID, file("${sampleID}.vcf") into vcf_ch

    shell:
	'''
    samtools faidx !{genome}
    bcftools mpileup -R !{bed} -f !{genome} !{bam} | bcftools call -mv -o !{sampleID}.vcf
    for sample in `bcftools query -l !{sampleID}.vcf`; do
        bcftools view -c1 -Oz -s $sample -o $sample.vcf !{sampleID}.vcf
    done
    '''
}


process NCM_run {

    cpus 1
    memory '16 GB'
    time { (2.hour + (2.hour * task.attempt)) }

    errorStrategy 'retry'
    maxRetries 3

    publishDir "$output", mode: 'copy'
    
    input:
    file (vcf) from vcf_ch.collect()

	output:
    file ("NCM_output") into ncm_ch
	
    script:
	"""
	ls \$PWD/*.vcf > listVCF

	mkdir NCM_output

	python \$NCM_HOME/ncm.py -V -l listVCF -bed ${bed} -O ./NCM_output
	Rscript NCM_output/r_script.r
    Rscript !{baseDir}/bin/plots.R $NCM_HOME 
	"""
}


process NCM_graphs {
    cpus 1
    memory '2 GB'

    publishDir "$output/NCM_output", mode: 'copy'
    
    
    when:
    params.NCM_labelfile!=null

    input:
    file ("NCM_output") from ncm_ch
    file labs from labelfile

	output:
    file ("*.xgmml") into graphs
	
    script:
	"""
    Rscript !{baseDir}/bin/plots.R $NCM_HOME !{labs}
	"""
}


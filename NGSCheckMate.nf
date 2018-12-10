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
    log.info "--input                  BAM FILES             BAM files (between quotes)"
    log.info "--output_dir             OUTPUT FOLDER         Output for NCM results"
    log.info "--ref_fasta              FASTA FILE            Reference FASTA file"
    log.info "--bed                    BED FILE              Selected SNPs file"
    exit 1
}


//
// Parameters Init
//
params.input         = null
params.output_dir    = "."
params.ref_fasta     = null
params.bed           = null


//
// Parse Input Parameters
//
bam_ch   = Channel
			.fromPath(params.input)
			.map { input -> tuple(input.baseName, input, input.parent / input.baseName + '.bai') }
output    = file(params.output_dir)
ref       = file(params.ref_fasta)
bed       = file(params.bed)


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

    output:
        set sampleID, file("${sampleID}.vcf") into vcf_ch

    script:
    """
        samtools faidx ${genome}
        
        bcftools mpileup -R ${bed} -f ${genome} ${bam} | bcftools call -mv -o ${sampleID}.vcf
    """

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

	"""
}


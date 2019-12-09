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
params.mem           = 16
params.cpu           = 4

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
    log.info "--input                  BAM FILES          List of BAM files (between quotes)"
    log.info "--output_folder          FOLDER             Output for NCM results"
    log.info "--ref                    FASTA FILE         Reference FASTA file"
    log.info "--bed                    BED FILE           Selected SNPs file"
    log.info "--NCM_labelfile          TSV FILE           tab-separated values file with 3 columns (vcf name, individual ID, sample ID) for generating xgmml graph file"
    log.info "--mem                    INTEGER            Memory (in GB)"
    log.info "--cpu                    INTEGER            Number of threads for germline calling"
    exit 0
}else{
  log.info "input_folder=${params.input_folder}"
  log.info "input=${params.input}"
  log.info "ref=${params.ref}"
  log.info "output_folder=${params.output_folder}"
  log.info "bed=${params.bed}"
  log.info "NCM_labelfile=${params.NCM_labelfile}"
  log.info "mem=${params.mem}"
  log.info "cpu=${params.cpu}"
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
ref       = file(params.ref)
bed       = file(params.bed)
labelfile = file(params.NCM_labelfile)
ncm_graphfiles = Channel.fromPath("$baseDir/bin/graph/*")

//
// Process Calling on SNP regions
//
process BCFTOOLS_calling{
    tag "$sampleID"

    cpus params.cpu
    memory params.mem+'G'

    input:
    file genome from ref 
    set sampleID, file(bam), file(bai) from bam_ch
    file bed

    output:
    file("*.vcf") into vcf_ch
    file("*.vcf.gz*") into vcfgz_ch

    publishDir params.output_folder+"/vcfs/", mode: 'copy'

    shell:
    cpus_mpileup = params.cpu.intdiv(2)
    cpus_call = params.cpu.intdiv(2)
	'''
    samtools faidx !{genome}
    bcftools mpileup --threads !{cpus_mpileup} --max-depth 5000 -Ou -I -R !{bed} -f !{genome} !{bam} | bcftools call --threads !{cpus_call} -c -o !{sampleID}_all.vcf
    for sample in `bcftools query -l !{sampleID}_all.vcf`; do
        bcftools view -Ou -s $sample !{sampleID}_all.vcf | bcftools sort -Ou | bcftools norm -d none -O v -o $sample.vcf
        bcftools view -Oz -o $sample.vcf.gz $sample.vcf
        tabix -p vcf $sample.vcf.gz
    done
    rm !{sampleID}_all.vcf
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

	output:
    file ("NCM_output") into ncm_ch
	
    script:
	"""
	ls \$PWD/*.vcf > listVCF

	mkdir NCM_output

	python \$NCM_HOME/ncm.py -V -l listVCF -bed ${bed} -O ./NCM_output
	"""
}

process NCM_graphs {
    cpus 1
    memory '2 GB'

    publishDir "${params.output_folder}/NCM_output", mode: 'copy'
    
    when:
    params.NCM_labelfile!=null

    input:
    file ("NCM_output") from ncm_ch
    file labs from labelfile
    file graphfiles from ncm_graphfiles.collect()

	output:
    file ("*.xgmml") into graphs
	
    when:
    labelfile.name != 'NO_FILE'

    shell:
	'''
    cat !{labs} | awk '{print $1".vcf\t"$2"\t"$3}' > input_plots.tsv
    Rscript !{baseDir}/bin/plots.R input_plots.tsv
	'''
}

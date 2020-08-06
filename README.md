# NGSCheckMate
## Nextflow pipeline to detect matched BAMs with [NGSCheckMate](https://github.com/parklab/NGSCheckMate).
[![CircleCI](https://circleci.com/gh/IARCbioinfo/NGSCheckMate-nf/tree/master.svg?style=svg)](https://circleci.com/gh/IARCbioinfo/NGSCheckMate-nf/tree/master)
[![Docker Hub](https://img.shields.io/badge/docker-ready-blue.svg)](https://hub.docker.com/repository/docker/iarcbioinfo/ngscheckmate-nf)
[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/4613)

<div style="text-align:center"><img src="https://camo.githubusercontent.com/371f23d984f8679c6562758f1e5b5e12397f1bef/68747470733a2f2f7061726b6c61622e6769746875622e696f2f4e4753436865636b4d6174652f6c6f676f2e737667" width="200" /></div>

## Description

Implementation of NGSCheckMate and its underlying subset calling, distibuted per sample.

## Dependencies 

1. Nextflow : for common installation procedures see the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository.
2. [NGSCheckMate](https://github.com/parklab/NGSCheckMate) (follow instructions, especially setting up `$NCM_HOME` variable)
3. [samtools](http://www.htslib.org/download/)
4. [bcftools](http://www.htslib.org/download/)

Additionally, the graph output option requires [R](https://cran.r-project.org/); see details below about this option.

## Input
  | Type      | Description     |
  |-----------|---------------|
  | --input   | your input BAM file(s) (do not forget the quotes e.g. `--input "test_*.bam"`). Warning : your BAM file(s) must be indexed, and the `test_*.bai` should be in the same folder.  |
  |  --input_folder  | Folder with BAM files  |
  | --input_file  | Input file (comma-separated) with 3 columns: ID (individual ID), suffix (suffix for sample names; e.g. RNA), and bam (path to bam file).|

A nextflow.config is also included, please modify it for suitability outside our pre-configured clusters ([see Nexflow configuration](https://www.nextflow.io/docs/latest/config.html#configuration-file)).

Note that the input_file format is tab-delimited text file; this file is used both to provide input bam file locations but also for the generation of the graphs.  The ID field must be unique to a subject (e.g. both tumor and normal samples from the same individual must have the same individual identifier). The bam field must be unique to a file name. For example, the following is a valid file:

ID  suffix  bam
NA06984 _RNA NA06984_T_transcriptome.bam\
NA06984 _WGS NA06984_T_genome.bam

## Parameters

  * #### Mandatory
| Name      | Example value | Description     |
|-----------|---------------|-----------------|
| --output_folder    |    results | the folder that will contain NGSCheckMate folder with all results in text files. |
| --ref    |          ref.fasta | your reference in FASTA |
| --bed |  SNP_GRCh38.bed | Panel of SNP bed file from [NGSCheckMate](https://github.com/parklab/NGSCheckMate/tree/master/SNP) |

Note that a bed file SNP_GRCh38.bed is provided, which is a liftOver of the files at https://github.com/parklab/NGSCheckMate/tree/master/SNP. To use other references, you can provide your own bedfile.


  * #### Optional
| Name      | Default value | Description     |
|-----------|---------------|-----------------|
| --mem   |   16 | Memory requested (in GB) for calling and NGSCheckmate run |
| --cpu    | 4 | Number of threads for germline calling |
|--bai_ext  | .bam.bai| Extenstion of bai files |

## Usage
  ```
  nextflow run NGSCheckMate-nf/ -r v1.1 -profile singularity --ref ref.fasta --input_folder BAM/
  ```
  
 To run the pipeline without singularity just remove "-profile singularity". Alternatively, one can run the pipeline using a docker container (-profile docker) the conda receipe containing all required dependencies (-profile conda).


## Output
  | Type      | Description     |
  |-----------|---------------|
  | vcfs    | a folder with the vcfs used for the matching |
  |  NCM_output/output*.txt   | NGSCheckmate output files with matches between files (see https://github.com/parklab/NGSCheckMate) |
  | NCM_output/output.pdf | hierarchical clustering plot from https://github.com/parklab/NGSCheckMate |
  | NCM_output/NCM_graph_wrongmatch.xgmml | graph with only the samples without a match (adapted from https://github.com/parklab/NGSCheckMate/blob/master/graph/ngscheckmate2xgmml.R) |
  | NCM_output/NCM_graph.xgmml | graph with all samples (adapted from https://github.com/parklab/NGSCheckMate/blob/master/graph/ngscheckmate2xgmml.R) |

Note that we recommend [Cytoscape](https://cytoscape.org/) to visualize the .xgmml graphs.

## Usage for Cobalt cluster
```
nextflow run iarcbioinfo/NGSCheckMate -profile cobalt --input "/data/test_*.bam" --output_dir /data/cohort_output --ref_fasta /ref/Homo_sapiens_assembly38.fasta --bed /home/user/bin/NGSCheckMate/SNP/SNP_GRCh38.bed
```

## FAQ

### Why are some files not included although the are in the intput_folder?
be careful that if bai files are missing for some bam files, the bam files will be ignored without the workflow returning an error

### What modifications have been done to the original NGSCheckMate code?
We provide a modified version of the graph/ngscheckmate2xgmml.R R script from https://github.com/parklab/NGSCheckMate to output graphs in .xgmml format. The modifications allow to represent all samples, even those that match, and improve a small glitch in the color palette.

## Contributions

  | Name      | Email | Description     |
  |-----------|---------------|-----------------| 
  | Nicolas Alcala*    | AlcalaN@iarc.fr    | Developer to contact for support |
  | Maxime Vallée |  | Developer |

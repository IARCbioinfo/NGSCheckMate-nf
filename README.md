# NGSCheckMate
Nextflow pipeline to detect matched BAMs with [NGSCheckMate](https://github.com/parklab/NGSCheckMate).

<div style="text-align:center"><img src="https://camo.githubusercontent.com/371f23d984f8679c6562758f1e5b5e12397f1bef/68747470733a2f2f7061726b6c61622e6769746875622e696f2f4e4753436865636b4d6174652f6c6f676f2e737667" width="200" /></div>

## Description

Implementation of NGSCheckMate and its underlying subset calling, distibuted per sample.

## Dependencies 

1. This pipeline is based on [nextflow](https://www.nextflow.io). As we have several nextflow pipelines, we have centralized the common information in the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository. Please read it carefully as it contains essential information for the installation, basic usage and configuration of nextflow and our pipelines.
2. [NGSCheckMate](https://github.com/parklab/NGSCheckMate) (follow instructions, especially setting up `$NCM_HOME` variable)
3. [samtools](http://www.htslib.org/download/)
4. [bcftools](http://www.htslib.org/download/)

## Input
  | Type      | Description     |
  |-----------|---------------|
  | --input   | your input BAM file(s) (do not forget the quotes e.g. `--input "test_*.bam"`). Warning : your BAM file(s) must be indexed, and the `test_*.bai` should be in the same folder.  |
  |  --input_folder  | Folder with BAM files  |

A nextflow.config is also included, please modify it for suitability outside our pre-configured clusters ([see Nexflow configuration](https://www.nextflow.io/docs/latest/config.html#configuration-file)).

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
| --NCM_labelfile | labels.tsv | tab-separated values file with 3 columns (vcf name, individual ID, sample ID) for generating xgmml graph file |

Note that the NCM_labelfile is optional; when provided, an extra step is computed using a modified version of the graph/ngscheckmate2xgmml.R R script from https://github.com/parklab/NGSCheckMate to output graphs in .xgmml format, which can be read by software [Cytoscape](https://cytoscape.org/). The format of file NCM_labelfile is similar to that of the original script from https://github.com/parklab/NGSCheckMate: a tab-delimited text file with 3 columns without header--a sample name (1st column) that must match the name in the SM field of the BAM header, an individual identifier (2nd column), and optionally, a file identifier (3rd column) for each line. The individual identifier must be unique to a subject (e.g. both tumor and normal samples from the same individual must have the same individual identifier). A file identifier must be unique to a file name. For example, the following is a valid file:

NA06984_T NA06984 NA06984_T_WES\
NA06984_N NA06984 NA06984_N_WES

where the bam files NA06984_T.bam and NA06984_N.bam respectively contain the following lines in their headers:\
@RG	ID:SRR098409	LB:Catch-36593	SM:NA06984_T	PI:110	CN:BI	PL:ILLUMINA	DS:SRP004078\
@RG	ID:SRR098409	LB:Catch-36593	SM:NA06984_N	PI:110	CN:BI	PL:ILLUMINA	DS:SRP004078


## Usage
  ```
  nextflow run NGSCheckMate-nf/ --ref ref.fasta --bed NGSCheckMate-nf/SNP_GRCh38.bed --input_folder BAM/ --NCM_labelfile NCM_labelfile.txt
  ```

## Output
  | Type      | Description     |
  |-----------|---------------|
  | vcfs    | a folder with the vcfs used for the matching |
  |  NCM_output/output*.txt   | NGSCheckmate output files with matches between files (see https://github.com/parklab/NGSCheckMate) |
  | NCM_output/output.pdf | hierarchical clustering plot from https://github.com/parklab/NGSCheckMate |
  | NCM_output/NCM_graph_wrongmatch.xgmml | graph with only the samples without a match (adapted from https://github.com/parklab/NGSCheckMate/blob/master/graph/ngscheckmate2xgmml.R) |
  | NCM_output/NCM_graph.xgmml | graph with all samples (adapted from https://github.com/parklab/NGSCheckMate/blob/master/graph/ngscheckmate2xgmml.R) |

Note that we recommend cytoscape to visualize the .xgmml graphs.


## Usage for Cobalt cluster
```
nextflow run iarcbioinfo/NGSCheckMate -profile cobalt --input "/data/test_*.bam" --output_dir /data/cohort_output --ref_fasta /ref/Homo_sapiens_assembly38.fasta --bed /home/user/bin/NGSCheckMate/SNP/SNP_GRCh38.bed
```

## FAQ

### Why are some files not included although the are in the intput_folder?
be careful that if bai files are missing for some bam files, the bam files will be ignored without the workflow returning an error

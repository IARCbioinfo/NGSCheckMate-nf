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

- `--input` : your input BAM file(s) (do not forget the quotes e.g. `--input "test_*.bam"`)
- `--output_dir` : the folder that will contain NGSCheckMate folder with all results in text files.
- `--ref_fasta` : your reference in FASTA. 
- `--bed` : Panel of SNP bed file from [NGSCheckMate](https://github.com/parklab/NGSCheckMate/tree/master/SNP). 


A nextflow.config is also included, please modify it for suitability outside our pre-configured clusters ([see Nexflow configuration](https://www.nextflow.io/docs/latest/config.html#configuration-file)).

## Usage for Cobalt cluster
```
nextflow run iarcbioinfo/NGSCheckMate.nf -profile cobalt --input "/data/test_*.bam" --output_dir /data/cohort_output --ref_fasta /ref/Homo_sapiens_assembly38.fasta --bed /home/user/bin/NGSCheckMate/SNP/SNP_GRCh38.bed
```


################## BASE IMAGE #####################
FROM continuumio/miniconda3:4.7.12

################## METADATA #######################

LABEL base_image="continuumio/miniconda3"
LABEL version="4.7.12"
LABEL software="NGSCheckMate-nf"
LABEL software.version="1.1"
LABEL about.summary="Container image containing all requirements for NGSCheckMate-nf"
LABEL about.home="http://github.com/IARCbioinfo/NGSCheckMate-nf"
LABEL about.documentation="http://github.com/IARCbioinfo/NGSCheckMate-nf/README.md"
LABEL about.license_file="http://github.com/IARCbioinfo/NGSCheckMate-nf/LICENSE.txt"
LABEL about.license="GNU-3.0"


################## MAINTAINER ######################
MAINTAINER **nalcala** <**alcalan@iarc.fr**>


################## INSTALLATION ######################
COPY environment.yml /
RUN apt-get update && apt-get install -y procps && apt-get clean -y
RUN conda env create -n ngscheckmate-nf -f /environment.yml && conda clean -a
RUN cd; git clone https://github.com/parklab/NGSCheckMate.git; echo "SAMTOOLS=samtools" >> NGSCheckMate/ncm.conf; echo "BCFTOOLS=bcftools" >> NGSCheckMate/ncm.conf
ENV NCM_HOME ~/NGSCheckMate
ENV PATH /opt/conda/envs/ngscheckmate-nf/bin:$PATH

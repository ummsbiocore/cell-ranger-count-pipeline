FROM ubuntu:22.04

LABEL description="Docker image containing all requirements via pipeline"

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

RUN apt-get update --fix-missing && \
    apt-get install -y vim wget bzip2 unzip ca-certificates curl git libtbb-dev gcc g++ libcairo2-dev pandoc libhdf5-dev cmake
    
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

RUN apt-get update && \
    curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64-2.0.30.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && ./aws/install && aws --version

COPY environment.yml /
RUN . /opt/conda/etc/profile.d/conda.sh && \ 
    conda activate base && \
    conda update conda && \
    conda env create -f /environment.yml && \
    conda clean -a

RUN mkdir -p /project /nl /mnt /share /pi
ENV PATH /opt/conda/envs/via/bin:$PATH
RUN git clone https://github.com/dolphinnext/tools /usr/local/bin/dolphin-tools
ENV PATH /opt/conda/envs/dolphinnext/bin:/usr/local/bin/dolphin-tools/:$PATH

## Cell Ranger ##
RUN cd /usr/bin && \
    wget https://web.dolphinnext.com/umw_biocore/dnext_data/bin/cellranger/cellranger-9.0.0.tar.gz && \
    tar -xzvf cellranger-9.0.0.tar.gz && \
    wget https://web.dolphinnext.com/umw_biocore/dnext_data/bin/cellranger/cellranger-atac-2.1.0.tar.gz && \
    tar -xzvf cellranger-atac-2.1.0.tar.gz
ENV PATH /usr/bin/cellranger-9.0.0:/usr/bin/cellranger-atac-2.1.0:$PATH

RUN pip install wheel
RUN pip install pandas
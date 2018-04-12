FROM ubuntu:trusty
MAINTAINER John Pellman <john.pellman@childmind.org>

ENV AFNIPATH /opt/afni
ENV PATH /code:/opt/afni:/usr/local/bin/miniconda/bin:${PATH}

# install dependencies
RUN apt-get update && \
    apt-get install -y pkg-config graphviz gsl-bin \
                       libexpat1-dev libgiftiio-dev libglu1-mesa libglu1-mesa-dev \
                       libgsl0-dev libjpeg-progs libxml2 libxml2-dev libxext-dev \
                       libxpm-dev libxp6 libxp-dev mesa-common-dev mesa-utils \
                       netpbm libpng-dev libfreetype6-dev libxml2-dev libxslt1-dev python-dev \
                       build-essential g++ libxft2 curl

RUN curl https://afni.nimh.nih.gov/pub/dist/tgz/linux_openmp_64.tgz -o /tmp/linux_openmp_64.tgz && \
    tar xzvf /tmp/linux_openmp_64.tgz -C /opt && mv /opt/linux_openmp_64 /opt/afni

# install miniconda
RUN curl -s https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -o /tmp/Miniconda2-latest-Linux-x86_64.sh && \
    bash /tmp/Miniconda2-latest-Linux-x86_64.sh -b -p /usr/local/bin/miniconda 
    
# install python requirements
RUN conda install -y pip scipy ipython
RUN pip install prov==1.5.0 networkx==1.11 nipype==0.13.1 python-dateutil==2.6.1 nibabel nitime pyyaml pandas seaborn pyPdf2 xhtml2pdf indi-tools configparser
    
COPY . /code

RUN cd /code && \
    pip install -e .
    # python setup.py build && python setup.py install

ENTRYPOINT [ "qap_run.py" ]

# docker run -ti --rm \
#     -v /tmp/workdir_ms \
#     -v /home/anibalsolon/examples:/dataset \
#     cmi/qap --data_config_file /dataset/bids_list.yml --pipeline_config_file /tmp/qap/configs/qap_config_ds003_ms.yml /dataset/ds003_downsampled /dataset/output participant
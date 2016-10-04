FROM ubuntu:trusty
MAINTAINER John Pellman <john.pellman@childmind.org>

ENV AFNIPATH /opt/afni/bin/
ENV PATH /code:/opt/afni/bin:/usr/local/bin/miniconda/bin:${PATH}

# install dependencies
RUN apt-get update && apt-get install -y wget
RUN apt-get install -y pkg-config graphviz gsl-bin \
    libexpat1-dev libgiftiio-dev libglu1-mesa libglu1-mesa-dev \
    libgsl0-dev libjpeg-progs libxml2 libxml2-dev libxext-dev \
    libxpm-dev libxp6 libxp-dev mesa-common-dev mesa-utils \
    netpbm libpng-dev libfreetype6-dev libxml2-dev libxslt1-dev python-dev \
    build-essential g++ libxft2

# install miniconda
RUN wget http://repo.continuum.io/miniconda/Miniconda-3.8.3-Linux-x86_64.sh && \
    bash Miniconda-3.8.3-Linux-x86_64.sh -b -p /usr/local/bin/miniconda && \
    rm -rf Miniconda-3.8.3-Linux-x86_64.sh && python -v

# install python requirements
RUN conda install -y pip scipy
RUN pip install nipype nibabel nitime pyyaml pandas seaborn pyPdf2 xhtml2pdf indi-tools

RUN wget http://afni.nimh.nih.gov/pub/dist/tgz/linux_openmp_64.tgz && \
    tar xzvf linux_openmp_64.tgz && mkdir -p /opt/afni && \
    mv linux_openmp_64/ /opt/afni/bin && \
    rm -rf linux_openmp_64.tgz

COPY . /tmp/qap

RUN cd /tmp/qap && \
    python setup.py build && python setup.py install

# don't forget ipython!
RUN conda install ipython

FROM afni/afni:latest

MAINTAINER Cameron Craddock <cameron.craddock@gmail.com>

ENV PATH /code:/usr/local/bin/miniconda/bin:${PATH}

COPY . /code

# install system packages
RUN apt-get update -y && \
    apt-get install -y xvfb

# install miniconda
RUN curl -s https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/Miniconda3-latest-Linux-x86_64.sh && \
    bash /tmp/Miniconda3-latest-Linux-x86_64.sh -b -p /usr/local/bin/miniconda

# install python requirements
RUN conda install -y pip scipy ipython

# install requirements and qap
RUN which pip && \
    cd /code && \
    pip install -r /code/requirements.txt && \
    pip install -e .

WORKDIR "/work"

ENTRYPOINT [ "qap_run.py" ]

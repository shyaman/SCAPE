FROM python:3.8.18-slim-bullseye

RUN apt-get update && apt-get install -y git wget build-essential zlib1g-dev libbz2-dev liblzma-dev r-base

# install bedtools v2.26.0
Run wget https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-2.26.0.tar.gz
RUN tar -zxvf bedtools-2.26.0.tar.gz
RUN cd bedtools2 && make && make install

# copy files to the docker image
COPY . /SCAPE

# install python dependencies
RUN cd SCAPE && pip install -r requirements.txt




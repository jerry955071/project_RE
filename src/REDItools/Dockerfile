# set base-image
FROM python:2.7.12

# select debian frontend
ARG DEBIAN_FRONTEND=noninteractive

# install dependencies
RUN pip install numpy pysam==0.15 && pip install fisher==0.1.4
# RUN apt-get update && apt-get install --force-yes samtools tabix

# install REDItools
RUN git clone https://github.com/BioinfoUNIBA/REDItools && \
    cd REDItools && \
    python setup.py install

# create mount point
RUN mkdir /data
WORKDIR /data

CMD REDItoolDnaRna.py



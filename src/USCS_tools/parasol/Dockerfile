FROM ubuntu:18.04

# Set working directory at root
WORKDIR /

# Install dependencies
RUN apt-get update \
  && apt-get install -y wget rsync yum git iproute2 libcurl4 libcurl4-openssl-dev\
  && rm -rf /var/lib/apt/lists/*

# Copy '/sbin/ip' to '/usr/sbin/ip' for /data/parasol/nodeInfo/nodeReport.sh
RUN cp -r /sbin/ip /usr/sbin/ip


# Install parasol job control system
RUN wget http://genomewiki.ucsc.edu/images/e/ea/ParasolInstall.sh.txt && \
    chmod +x /ParasolInstall.sh.txt && \
    /ParasolInstall.sh.txt

# Add parasol tool sets to path
ENV PATH="/data/bin:/data/scripts:${PATH}"

# Initialize node instances
# SSH key (not effecting)
# RUN echo | ssh-keygen -N "" -t rsa -q && \
#    cat $HOME/.ssh/id_rsa.pub >> $HOME/.ssh/authorized_keys
   
# usage: nodeReport.sh N nodeInfoDir
# use N=2 for a shared computer to reserve 2 CPUs for OS and the parasol hub process.
# RUN /data/parasol/nodeInfo/nodeReport.sh 2 /data/parasol/nodeInfo

# /data/parasol/initParasol initialize

RUN apt-get update && \
  apt-get install -y openssh-server && \
  rm -rf /var/lib/apt/lists/*

RUN echo export PATH="/data/bin:/data/scripts:${PATH}" >> /root/.bashrc


RUN apt-get update && \
  apt-get install -y csh && \
  rm -rf /var/lib/apt/lists/*

COPY myDoBlastzChainNet.pl /data/scripts
RUN chmod a+x /data/scripts/myDoBlastzChainNet.pl
# service ssh start && \
#   echo | ssh-keygen -N "" -t rsa -q && \
#   cat $HOME/.ssh/id_rsa.pub >> $HOME/.ssh/authorized_keys && \
#   /data/parasol/nodeInfo/nodeReport.sh 2 /data/parasol/nodeInfo && \
#   /data/parasol/initParasol initialize

COPY gensub2 /data/bin

COPY blastz-run-ucsc /data/scripts/

RUN chmod a+x /data/scripts/blastz-run-ucsc
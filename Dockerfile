FROM debian

MAINTAINER Chieh Lin <chiehl1@cs.cmu.edu>

RUN apt-get update
RUN apt-get -y install python2.7 python-pip
RUN apt-get -y install graphviz
RUN apt-get -y install vim
RUN apt-get -y install curl
RUN curl https://bootstrap.pypa.io/get-pip.py | python
RUN ln -s /usr/local/bin/pip /usr/bin/pip
RUN pip install --upgrade scdiff
RUN pip install cvxpy
RUN pip install progressbar
RUN apt-get -y install libgraphviz-dev 
RUN pip install pygraphviz
RUN pip install networkx
RUN pip install pandas
RUN pip install argparse
RUN pip install sklearn
RUN apt-get -y install r-base
RUN apt-get -y install libssl-dev
RUN apt-get -y install libcurl4-openssl-dev
RUN R --vanilla -e 'install.packages("devtools",repos="https://cran.cnr.berkeley.edu")'
RUN R --vanilla -e 'library(devtools);install_github("statsmaths/genlasso")'
RUN pip install rpy2==2.8.6
RUN apt-get -y install git
RUN apt-get -y install python-tk
RUN pip install -U statsmodels

RUN git clone https://github.com/jessica1338/CSHMM-TF-for-time-series-scRNA-Seq.git
RUN mkdir CSHMM-TF-for-time-series-scRNA-Seq/TF_analysis_result

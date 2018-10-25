FROM debian

MAINTAINER Chieh Lin <chiehl1@cs.cmu.edu>

RUN apt-get update
RUN apt-get install python2.7
RUN apt-get install graphviz
RUN apt-get install vim
RUN apt-get install curl
RUN curl https://bootstrap.pypa.io/get-pip.py | python
RUN ln -s /usr/local/bin/pip /usr/bin/pip
RUN pip install --upgrade scdiff
RUN pip install cvxpy
RUN pip install progressbar
RUN apt-get instaall libgraphviz-dev 
RUN apt-get install libgraphviz-dev 
RUN pip install pygraphviz
RUN pip install networkx
RUN pip install pandas
RUN pip install argparse
RUN pip install sklearn
RUN apt-get install r-base
RUN apt-get install libssl-dev
RUN apt-get install libcurl4-openssl-dev
RUN R --vanilla -e 'install.packages("devtools",repos="https://cran.cnr.berkeley.edu")'
RUN R --vanilla -e 'library(devtools);install_github("statsmaths/genlasso")'
RUN pip install rpy2==2.8.6



## install scdiff 

## install genlasso
## install graphviz
## install python 2.7
## install cvxpy
## install progressbar
## install pygraphviz
## install matplotlib
## install networkx
## install numpy
## install pandas
## install argparse
## install sklearn



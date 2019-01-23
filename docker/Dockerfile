FROM python:3.6-stretch

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install apt-utils -y 

RUN apt-get install -y net-tools \
                       build-essential \
                       cmake && \
    mkdir -p /gurobi /opt/gurobi

COPY gurobi7.5.2_linux64.tar.gz /gurobi/gurobi7.5.2_linux64.tar.gz

RUN tar xvf /gurobi/gurobi7.5.2_linux64.tar.gz -C /gurobi && \
    cd /gurobi/gurobi752/linux64/lib && \
	rm libgurobi_c++.a && \
	ln -s libgurobi_g++5.2.a libgurobi_c++.a 

RUN echo "TOKENSERVER=license" > /opt/gurobi/gurobi.lic
    
RUN cd tmp && \
    wget http://lemon.cs.elte.hu/pub/sources/lemon-1.3.1.tar.gz && \
    tar xvf lemon-1.3.1.tar.gz && cd lemon-1.3.1 && \
    mkdir build && cd build && mkdir /lemon && \
    cmake CMAKE_INSTALL_PREFIX=/lemon .. && make && make install && \ 
    cd / && rm -r tmp/lemon-1.3.1*

RUN git clone https://github.com/sebwink/deregnet

COPY common.mak common.mak
COPY gurobi_version.mak gurobi_version.mak

RUN cp common.mak /deregnet/common.mak && \
    cp gurobi_version.mak /deregnet/gurobi_version.mak && \
	cd /deregnet && \
	make && \
	ln -s /deregnet/bin/avgdrgnt /usr/local/bin/avgdrgnt && \
	ln -s /deregnet/bin/drgnt /usr/local/bin/drgnt && \
	ln -s /deregnet/bin/avgdrgnt.py /usr/local/bin/avgdrgnt.py && \
	ln -s /deregnet/bin/drgnt.py /usr/local/bin/drgnt.py && \
	ln -s /deregnet/python/deregnet /usr/local/lib/python3.6/site-packages/deregnet

RUN	python3 -m pip install --upgrade pip && \
	python3 -m pip install python-igraph && \
    python3 -m pip install pymongo && \
	python3 -m pip install pandas && \
	python3 -m pip install rarfile && \
	python3 -m pip install pyyaml && \
	python3 -m pip install requests

RUN mkdir -p /root/.deregnet/graphs/kegg && \
    mkdir /root/.deregnet/graphs/omnipath

COPY graphs/kegg/kegggraph/hsa/kegg_hsa.sif /root/.deregnet/graphs/kegg/kegg_hsa.sif
COPY graphs/kegg/kegggraph/mmu/kegg_mmu.sif /root/.deregnet/graphs/kegg/kegg_mmu.sif
COPY graphs/kegg/kegggraph/sce/kegg_sce.sif /root/.deregnet/graphs/kegg/kegg_sce.sif
COPY graphs/omnipath /root/.deregnet/graphs/omnipath

CMD /bin/bash

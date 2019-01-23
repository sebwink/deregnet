FROM sebwink/gurobi

RUN apt-get install -y build-essential cmake

RUN python3 -m pip install --upgrade pip && \
    python3 -m pip install python-igraph

RUN python3 -m pip install pandas

RUN cd tmp && \
    wget http://lemon.cs.elte.hu/pub/sources/lemon-1.3.1.tar.gz && \
    tar xvf lemon-1.3.1.tar.gz && cd lemon-1.3.1 && \
    mkdir build && cd build && mkdir -p /opt/lemon/1.3.1 && \
    cmake -DCMAKE_INSTALL_PREFIX=/opt/lemon/1.3.1 .. && make && make install && \ 
    cd / && rm -r tmp/lemon-1.3.1*

RUN mkdir /deregnet
WORKDIR /deregnet 

COPY common.mak .
COPY gurobi_version.mak .
COPY drgnt.mak .
COPY avgdrgnt.mak .
COPY build.sh .
COPY Makefile .
COPY include include
COPY src src
COPY grbfrc grbfrc
COPY python python
COPY bin bin

RUN cd /deregnet && \
    make && \
    ln -s /deregnet/bin/avgdrgnt /usr/local/bin/avgdrgnt && \
    ln -s /deregnet/bin/drgnt /usr/local/bin/drgnt && \
    ln -s /deregnet/bin/avgdrgnt.py /usr/local/bin/avgdrgnt.py && \
    ln -s /deregnet/bin/drgnt.py /usr/local/bin/drgnt.py && \
    ln -s /deregnet/python/deregnet /usr/local/lib/python3.6/site-packages/deregnet

RUN python3 -m pip install biomap-utils

RUN mkdir /io
WORKDIR /io

ENTRYPOINT ["/deregnet/bin/avgdrgnt.py"]

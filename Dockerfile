ARG GRBFRC_IMAGE_TAG
ARG GRBFRC_GUROBI_VERSION

FROM sebwink/lemon-headers:131 as lemon

FROM sebwink/libgrbfrc-grb${GRBFRC_GUROBI_VERSION}:${GRBFRC_IMAGE_TAG}

ARG GUROBI_USER
ARG GRBFRC_GUROBI_VERSION

# CMAKE version
ARG CMAKE_VERSION
ENV CMAKE_VERSION=${CMAKE_VERSION:-3.20.1}

COPY --from=lemon /usr/local/include/lemon /usr/local/include/lemon

USER root

RUN apt-get update && \
    apt-get upgrade -y && \
	apt-get install -y build-essential && \
	apt-get install -y python3-dev python3-pip && \
	mkdir /deregnet
	
# CMAKE
RUN wget https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}-linux-x86_64.tar.gz && \
    tar xvf cmake-${CMAKE_VERSION}-linux-x86_64.tar.gz && rm cmake-${CMAKE_VERSION}-linux-x86_64.tar.gz && \
    cd cmake-${CMAKE_VERSION}-linux-x86_64 && \
    cp bin/* /usr/local/bin/ && \
    mkdir -p /usr/local/man/man1 && \
    cp man/man1/* /usr/local/man/man1/ && \
    mkdir -p /usr/local/man/man7 && \
    cp man/man7/* /usr/local/man/man7/ && \
    mkdir -p /usr/local/doc && \
    cp -r doc/cmake /usr/local/doc/ && \
    mkdir -p /usr/local/share/bash-completion/completions && \
    cp share/bash-completion/completions/* /usr/local/share/bash-completion/completions/ && \
    CMAKE_MAJOR_MINOR=$(echo $CMAKE_VERSION | sed "s/\([0-9]\+\.[0-9]\+\)\.[0-9]\+/\1/g") && \
    cp -r share/cmake-${CMAKE_MAJOR_MINOR} /usr/local/share/cmake-${CMAKE_MAJOR_MINOR} && \
    cd .. && rm -r cmake-${CMAKE_VERSION}-linux-x86_64

WORKDIR /deregnet 

COPY Makefile .
COPY src src
COPY upstream/libgrbfrc/gurobi.mak gurobi.mak

RUN GUROBI_VERSION_TAG=${GRBFRC_GUROBI_VERSION} LIBGRBFRC=/usr/local/include GUROBI_MAKEFILE=gurobi.mak make all

RUN python3 -m pip install pandas && \
	python3 -m pip install biomap-utils && \
	apt-get install -y libz-dev && \
	apt-get install -y libxml2-dev && \
	apt-get install -y git && \
	apt-get install -y libtool && \
	apt-get install bison -y && \
        apt-get install byacc -y && \
	apt-get install flex -y && \
        python3 -m pip install python-igraph && \
	python3 -m pip install ipython jupyterlab

COPY python python 

WORKDIR /io

RUN chown -R ${GUROBI_USER}:${GUROBI_USER} /io && \
		echo "GUROBI_USER: ${GUROBI_USER}" && \
		mkdir -p /home/${GUROBI_USER} && \
		chown -R ${GUROBI_USER}:${GUROBI_USER} /home/${GUROBI_USER}

RUN ln -s /deregnet/python/deregnet /usr/local/lib/python3.7/dist-packages/deregnet

USER ${GUROBI_USER}

ENV PATH=/deregnet/bin:${PATH}
ENV PATH=/deregnet/python/scripts:${PATH}

ENTRYPOINT ["sh", "-c"]

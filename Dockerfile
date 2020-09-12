ARG GRBFRC_IMAGE_TAG=latest
ARG GRBFRC_GUROBI_VERSION=9.0.2

FROM sebwink/lemon-headers:131 as lemon

FROM sebwink/libgrbfrc-grb${GRBFRC_GUROBI_VERSION}:${GRBFRC_IMAGE_TAG}

ARG GUROBI_USER

COPY --from=lemon /usr/local/include/lemon /usr/local/include/lemon

USER root

RUN apt-get update && \
    apt-get upgrade -y && \
	apt-get install -y build-essential && \
	apt-get install -y python3-dev python3-pip && \
	mkdir /deregnet

WORKDIR /deregnet 

COPY Makefile .
COPY src src
COPY upstream/libgrbfrc/gurobi.mak gurobi.mak

RUN LIBGRBFRC=/usr/local/include GUROBI_MAKEFILE=gurobi.mak make all

RUN python3 -m pip install pandas && \
	python3 -m pip install biomap-utils && \
	apt-get install -y libz-dev && \
	apt-get install -y libxml2-dev && \
	apt-get install -y git && \
	apt-get install -y libtool && \
	apt-get install bison -y && \
    apt-get install byacc -y && \
	apt-get install flex -y && \
    python3 -m pip install python-igraph

COPY python python 

WORKDIR /io

RUN chown -R ${GUROBI_USER}:${GUROBI_USER} /io && \
		echo "GUROBI_USER: ${GUROBI_USER}" && \
		mkdir -p /home/${GUROBI_USER} && \
		chown -R ${GUROBI_USER}:${GUROBI_USER} /home/${GUROBI_USER}

USER ${GUROBI_USER}

ENV PATH=/deregnet/python/scripts:${PATH}

ENTRYPOINT ["/deregnet/python/scripts/avgdrgnt.py"]

FROM debian:stretch-slim as setup

ARG GUROBI_MAJOR_VERSION=8.1
ARG GUROBI_VERSION=8.1.0
ARG GUROBI_TAG=810

RUN apt-get update && \
    apt-get upgrade -y && \
	apt-get install -y wget && \
	wget http://packages.gurobi.com/${GUROBI_MAJOR_VERSION}/gurobi${GUROBI_VERSION}_linux64.tar.gz && \
    tar xvfz gurobi${GUROBI_VERSION}_linux64.tar.gz
    

FROM debian:stretch-slim

ARG GUROBI_VERSION=8.1.0
ARG GUROBI_TAG=810

COPY --from=setup /gurobi${GUROBI_TAG}/linux64/lib/libgurobi.so.${GUROBI_VERSION} /usr/local/lib/
COPY --from=setup /gurobi${GUROBI_TAG}/linux64/lib/libgurobi_g++5.2.a /usr/local/lib/libgurobi_c++.a
COPY --from=setup /gurobi${GUROBI_TAG}/linux64/include/gurobi_c.h /usr/local/include/
COPY --from=setup /gurobi${GUROBI_TAG}/linux64/include/gurobi_c++.h /usr/local/include/

RUN mkdir /opt/gurobi && \
    echo "TOKENSERVER=license" > /opt/gurobi/gurobi.lic && \
	ln -s /usr/local/lib/libgurobi.so.${GUROBI_VERSION} /usr/local/lib/libgurobi${GUROBI_TAG}.so && \
    ln -s /usr/local/lib/libgurobi${GUROBI_TAG}.so /usr/local/lib/libgurobi.so

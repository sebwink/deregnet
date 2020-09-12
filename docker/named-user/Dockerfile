ARG BASEIMAGE

FROM ${BASEIMAGE}

ARG GUROBI_USER

USER root 

RUN groupadd ${GUROBI_USER} && useradd -r -m -g ${GUROBI_USER} ${GUROBI_USER} && \
		chown -R ${GUROBI_USER}:${GUROBI_USER} /deregnet && \
		chown -R ${GUROBI_USER}:${GUROBI_USER} /io

USER ${GUROBI_USER}

WORKDIR /io

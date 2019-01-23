FROM python:3.6-stretch

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install apt-utils -y 

ARG GUROBI_MAJOR_VERSION=8.1
ARG GUROBI_VERSION=8.1.0
ARG GUROBI_TAG=810

RUN mkdir -p /opt/gurobi/io

RUN wget http://packages.gurobi.com/${GUROBI_MAJOR_VERSION}/gurobi${GUROBI_VERSION}_linux64.tar.gz && \
    tar xvfz gurobi${GUROBI_VERSION}_linux64.tar.gz -C /opt/gurobi && \
    rm -r gurobi${GUROBI_VERSION}_linux64.tar.gz

RUN mv /opt/gurobi/gurobi${GUROBI_TAG}/linux64/* /opt/gurobi && \
    rm -r /opt/gurobi/gurobi${GUROBI_TAG} && \
    cd /opt/gurobi/lib && \
    rm libgurobi_c++.a && \
    ln -s libgurobi_g++5.2.a libgurobi_c++.a 

RUN echo "TOKENSERVER=license" > /opt/gurobi/gurobi.lic

RUN python3 -m pip install --upgrade pip && \
    cd /opt/gurobi && \
    python3 setup.py install

ENV GUROBI_HOME=/opt/gurobi
ENV LD_LIBRARY_PATH=/opt/gurobi/lib:$LD_LIBRARY_PATH
ENV PATH=/opt/gurobi/bin:$PATH

COPY cmd.sh .

CMD ["bash", "cmd.sh"]

FROM sebwink/lemon-headers:131 as lemon

FROM sebwink/grbfrc as build

COPY --from=lemon /usr/local/include/lemon /usr/local/include/lemon

RUN apt-get update && \
    apt-get upgrade -y && \
	apt-get install -y build-essential && \
    mkdir -p /deregnet/bin /deregnet/build
WORKDIR /deregnet 

COPY docker/deregnet/common.mak .
COPY docker/deregnet/drgnt.mak .
COPY docker/deregnet/avgdrgnt.mak .
COPY include include
COPY src src

RUN make -f avgdrgnt.mak && \
    make -f drgnt.mak

RUN python3 -m pip install --upgrade pip && \
	python3 -m pip install pandas && \
	python3 -m pip install biomap-utils && \
	apt-get install -y libz-dev && \
	apt-get install -y libxml2-dev && \
    python3 -m pip install python-igraph

FROM python:3.6-slim

RUN apt-get update && \
    apt-get upgrade -y && \
	apt-get install -y libxml2 libz-dev

COPY --from=build /usr/local/lib /usr/local/lib
COPY --from=build /deregnet/bin/drgnt /usr/local/bin/
COPY --from=build /deregnet/bin/avgdrgnt /usr/local/bin/
COPY --from=build /usr/local/lib/python3.6/site-packages /usr/local/lib/python3.5/site-packages

COPY python/deregnet /usr/local/lib/python3.6/site-packages/deregnet 
COPY bin/drgnt.py /usr/local/bin/
COPY bin/avgdrgnt.py /usr/local/bin/

WORKDIR /io

ENTRYPOINT ["avgdrgnt.py"]

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

RUN apt-get install -y python3 && \
	apt-get install -y python3-pip && \
	python3 -m pip install --upgrade pip && \
	python3 -m pip install pandas && \
	python3 -m pip install biomap-utils && \
	apt-get install -y libz-dev && \
	apt-get install -y libxml2-dev && \
    python3 -m pip install python-igraph

FROM debian:stretch-slim

RUN apt-get update && \
    apt-get upgrade -y && \
	apt-get install -y python3-minimal libxml2 libz-dev

COPY --from=build /usr/local/lib /usr/local/lib
COPY --from=build /deregnet/bin/drgnt /usr/local/bin/
COPY --from=build /deregnet/bin/avgdrgnt /usr/local/bin/
COPY --from=build /usr/local/lib/python3.5/dist-packages /usr/local/lib/python3.5/dist-packages

COPY python/deregnet /usr/local/lib/python3.5/dist-packages/deregnet 
COPY bin/drgnt.py /usr/local/bin/
COPY bin/avgdrgnt.py /usr/local/bin/

RUN python3 -m pip install six && mkdir /io
WORKDIR /io

ENTRYPOINT ["avgdrgnt.py"]

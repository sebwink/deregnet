FROM debian:buster-slim as build

RUN apt-get update && \
    apt-get upgrade -y && \
	apt-get install -y wget cmake build-essential && \
    wget http://lemon.cs.elte.hu/pub/sources/lemon-1.3.1.tar.gz && \
    tar xvf lemon-1.3.1.tar.gz && \
	cd lemon-1.3.1 && \
    mkdir build && cd build && \
    cmake .. && make && make install

FROM debian:stretch-slim 

COPY --from=build /usr/local/include /usr/local/include

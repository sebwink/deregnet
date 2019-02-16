FROM sebwink/gurobi-cpp:810 as build 

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y build-essential && \
    mkdir grbfrc   

WORKDIR /grbfrc

COPY docker/Makefile .
COPY src src
COPY include include

RUN make 

FROM sebwink/gurobi-cpp:810

COPY --from=build /grbfrc/lib/libgrbfrc.so /usr/local/lib/
COPY --from=build /grbfrc/include/grbfrc /usr/local/include/grbfrc

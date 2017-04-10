.PHONY: destroy 

all: lib/libderegnet.so

include common.mak

OBJECTS=build/AvgStartHeuristic.o \
        build/AvgSuboptimalStartHeuristic.o \
		build/DeregnetData.o \
		build/DeregnetStartHeuristic.o \
		build/LazyConstraintCallback.o \
		build/StartHeuristic.o \
		build/SuboptimalStartHeuristic.o \
		build/utils.o

LDFLAGS=-L${GUROBI_HOME}/lib -lgurobi_c++ -lgurobi${GUROBI_VERSION_SUFFIX} -L${GRBFRC_HOME}/lib -lgrbfrc
RUNPATH=${GUROBI_HOME}/lib:../${GRBFRC_HOME}/lib

lib/libderegnet.so : 
	$(CXX) -shared -o $@ $(OBJECTS) $(LDFLAGS) -Wl,-rpath=$(RUNPATH)
	
destroy:
	rm -f lib/*

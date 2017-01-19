CXX=c++
CXXFLAGS=-std=c++11 -Wall -ggdb3

LEMON_VERSION=1.3.1
LEMON_HOME=/opt/lemon/${LEMON_VERSION}

GUROBI_VERSION_MAJOR=7
GUROBI_VERSION_MINOR=0
GUROBI_VERSION=${GUROBI_VERSION_MAJOR}.${GUROBI_VERSION_MINOR}.1
GUROBI_VERSION_SUFFIX=${GUROBI_VERSION_MAJOR}${GUROBI_VERSION_MINOR}

GUROBI_HOME=/opt/gurobi/${GUROBI_VERSION}

INCLUDE=-Iinclude -I${LEMON_HOME}/include
LDFLAGS=-L${LEMON_HOME}/lib -lemon

build/utils.o : src/utils.cpp include/deregnet/utils.h include/deregnet/usinglemon.h
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(INCLUDE)

build/DeregnetStartHeuristic.o : src/DeregnetStartHeuristic.cpp include/deregnet/DeregnetStartHeuristic.h include/deregnet/usinglemon.h
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(INCLUDE)

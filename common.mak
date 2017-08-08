CXX=c++
CXXFLAGS=-std=c++11 -Wall -fPIC 

LEMON_VERSION=1.3.1
LEMON_HOME=/share/opt/x86_64_sl7/lemon

include gurobi_version.mak

GRBFRC_HOME=/home/winkler/prjcts/deregnet/grbfrc

INCLUDE=-Iinclude -I${LEMON_HOME}/include -I${GUROBI_HOME}/include -I${GRBFRC_HOME}/include

build/utils.o : src/utils.cpp include/deregnet/utils.h include/deregnet/usinglemon.h
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(INCLUDE)

build/DeregnetStartHeuristic.o : src/DeregnetStartHeuristic.cpp include/deregnet/DeregnetStartHeuristic.h include/deregnet/usinglemon.h
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(INCLUDE)

build/DeregnetData.o : src/DeregnetData.cpp include/deregnet/DeregnetData.h
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(INCLUDE)


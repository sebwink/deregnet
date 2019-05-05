CXX=c++
CXXFLAGS=-std=c++11 -Wall -fPIC 

INCLUDE=-Iinclude 

build/utils.o : src/utils.cpp include/deregnet/utils.h include/deregnet/usinglemon.h
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(INCLUDE)

build/DeregnetStartHeuristic.o : src/DeregnetStartHeuristic.cpp include/deregnet/DeregnetStartHeuristic.h include/deregnet/usinglemon.h
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(INCLUDE)

build/DeregnetData.o : src/DeregnetData.cpp include/deregnet/DeregnetData.h
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(INCLUDE)


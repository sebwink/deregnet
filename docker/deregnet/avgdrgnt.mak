.PHONY: all clean destroy

all : bin/avgdrgnt

include common.mak

LDFLAGS = -lgrbfrc -lgurobi_c++ -lgurobi

objects = build/utils.o \
          build/DeregnetData.o \
          build/DeregnetStartHeuristic.o \
          build/avgdrgnt.o \
          build/AvgStartHeuristic.o \
	      build/AvgSuboptimalStartHeuristic.o \
	      build/GrbfrcLazyConstraintCallback.o \
	      build/LazyConstraintCallback.o

bin/avgdrgnt : $(objects)
	$(CXX) -o $@ $^ $(LDFLAGS)

AVGDRGNT_INCLUDE_DEPS=include/deregnet/version.h \
                      include/deregnet/utils.h \
		              include/deregnet/DeregnetFinder.h \
		              include/deregnet/DeregnetModel.h \
		              include/deregnet/AvgdrgntData.h

build/avgdrgnt.o : src/bin/avgdrgnt.cpp $(AVGDRGNT_INCLUDE_DEPS)
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(INCLUDE)

# LazyConstraintCallback

GRBFRCLAZYCONSTRAINTCALLBACK_INCLUDE_DEPS=include/deregnet/usinglemon.h \
                                          include/deregnet/GrbfrcLazyConstraintCallback.h

build/GrbfrcLazyConstraintCallback.o : src/GrbfrcLazyConstraintCallback.cpp $(GRBFRCLAZYCONSTRAINTCALLBACK_INCLUDE_DEPS)
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(INCLUDE)

# LazyConstraintCallback

LAZYCONSTRAINTCALLBACK_INCLUDE_DEPS=include/deregnet/usinglemon.h \
                                    include/deregnet/LazyConstraintCallback.h

build/LazyConstraintCallback.o : src/LazyConstraintCallback.cpp $(LAZYCONSTRAINTCALLBACK_INCLUDE_DEPS)
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(INCLUDE)

# AvgStartHeuristic

AVGSTARTHEURISTIC_INCLUDE_DEPS=include/deregnet/DeregnetStartHeuristic.h \
                               include/deregnet/AvgStartHeuristic.h \
                               include/deregnet/usinglemon.h

build/AvgStartHeuristic.o : src/AvgStartHeuristic.cpp $(AVGSTARTHEURISTIC_INCLUDE_DEPS)
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(INCLUDE)

# SuboptimalStartHeuristic

AVGSUBOPTIMALSTARTHEURISTIC_INCLUDE_DEPS=include/deregnet/DeregnetStartHeuristic.h \
                                         include/deregnet/AvgSuboptimalStartHeuristic.h \
                                         include/deregnet/usinglemon.h

build/AvgSuboptimalStartHeuristic.o : src/AvgSuboptimalStartHeuristic.cpp $(AVGSUBOPTIMALSTARTHEURISTIC_INCLUDE_DEPS)
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(INCLUDE)


clean : 
	rm -f $(objects)

destroy :
	rm -f $(objects)
	rm -f bin/avgdrgnt

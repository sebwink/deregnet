.PHONY: all clean destroy

all : bin/drgnt

include common.mak

CXXFLAGS += -DRGNT_DEBUG

INCLUDE += -I${GUROBI_HOME}/include
LDFLAGS += -L${GUROBI_HOME}/lib -lgurobi_c++ -lgurobi${GUROBI_VERSION_SUFFIX}
RUNPATH = ${LEMON_HOME}/lib:${GUROBI_HOME}/lib

objects = build/utils.o \
          build/DeregnetData.o \
          build/LazyConstraintCallback.o \
		  build/DeregnetStartHeuristic.o \
		  build/StartHeuristic.o \
		  build/SuboptimalStartHeuristic.o \
		  build/drgnt.o

bin/drgnt : $(objects)
	$(CXX) -o $@ $^ $(LDFLAGS) -Wl,-rpath=$(RUNPATH)

DRGNT_INCLUDE_DEPS=include/deregnet/version.h \
                   include/deregnet/utils.h \
				   include/deregnet/DeregnetFinder.h \
				   include/deregnet/DrgntData.h

build/drgnt.o : src/bin/drgnt.cpp $(DRGNT_INCLUDE_DEPS)
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(INCLUDE)

# LazyConstraintCallback

LAZYCONSTRAINTCALLBACK_INCLUDE_DEPS=include/deregnet/usinglemon.h \
                                    include/deregnet/LazyConstraintCallback.h

build/LazyConstraintCallback.o : src/LazyConstraintCallback.cpp $(LAZYCONSTRAINTCALLBACK_INCLUDE_DEPS)
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(INCLUDE)

# StartHeuristic

STARTHEURISTIC_INCLUDE_DEPS=include/deregnet/DeregnetStartHeuristic.h \
                            include/deregnet/StartHeuristic.h \
							include/deregnet/usinglemon.h

build/StartHeuristic.o : src/StartHeuristic.cpp $(STARTHEURISTIC_INCLUDE_DEPS)
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(INCLUDE)

# SuboptimalStartHeuristic

SUBOPTIMALSTARTHEURISTIC_INCLUDE_DEPS=include/deregnet/DeregnetStartHeuristic.h \
                                      include/deregnet/SuboptimalStartHeuristic.h \
					                  include/deregnet/usinglemon.h

build/SuboptimalStartHeuristic.o : src/SuboptimalStartHeuristic.cpp $(SUBOPTIMALSTARTHEURISTIC_INCLUDE_DEPS)
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(INCLUDE)


clean : 
	rm -f $(objects)

destroy :
	rm -f $(objects)
	rm -f bin/drgnt

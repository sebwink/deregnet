.PHONY: all clean destroy

all : bin/drgnt

include common.mak

INCLUDE += -I${GUROBI_HOME}/include
LDFLAGS += -L${GUROBI_HOME}/lib -lgurobi_c++ -lgurobi${GUROBI_VERSION_SUFFIX}
RUNPATH = ${LEMON_HOME}/lib:${GUROBI_HOME}/lib

objects = build/utils.o build/SbgrphFinder.o build/SbgrphModel.o build/drgnt.o

bin/drgnt : $(objects)
	$(CXX) -o $@ $^ $(LDFLAGS) -Wl,-rpath=$(RUNPATH)

DRGNT_INCLUDE_DEPS=include/deregnet/version.h \
                   include/deregnet/utils.h \
				   include/deregnet/SbgrphFinder.h

build/drgnt.o : src/bin/drgnt.cpp $(DRGNT_INCLUDE_DEPS)
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(INCLUDE)

# SbgrphFinder

SBGRPHFINDER_INCLUDE_DEPS=include/deregnet/SbgrphFinder.h \
                          include/deregnet/SbgrphModel.h \
						  include/deregnet/utils.h \
						  include/deregnet/usinglemon.h

build/SbgrphFinder.o : src/SbgrphFinder.cpp $(SBGRPHFINDER_INCLUDE_DEPS)
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(INCLUDE)

# SbgrphModel

SBGRPHMODEL_INCLUDE_DEPS=$(SBGRPHFINDER_INCLUDE_DEPS)

build/SbgrphModel.o : src/SbgrphModel.cpp $(SBGRPHMODEL_INCLUDE_DEPS)
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(INCLUDE)


clean : 
	rm -f $(objects)

destroy :
	rm -f $(objects)
	rm -f bin/drgnt

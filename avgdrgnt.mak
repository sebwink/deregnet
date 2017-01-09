.PHONY: all clean destroy

all : bin/avgdrgnt

include common.mak

GRBFRC_HOME=/home/sebastian/wrk/grbfrc/grbfrc

INCLUDE += -I${GRBFRC_HOME}/include
LDFLAGS += -L${GRBFRC_HOME}/lib -lgrbfrc
RUNPATH = ${LEMON_HOME}/lib:${GRBFRC_HOME}/lib

objects = build/utils.o build/AvgSbgrphFinder.o build/AvgSbgrphModel.o build/avgdrgnt.o

bin/avgdrgnt : $(objects)
	$(CXX) -o $@ $^ $(LDFLAGS) -Wl,-rpath=$(RUNPATH)

DRGNT_INCLUDE_DEPS=include/deregnet/version.h \
                   include/deregnet/utils.h \
				   include/deregnet/AvgSbgrphFinder.h

build/avgdrgnt.o : src/bin/avgdrgnt.cpp $(DRGNT_INCLUDE_DEPS)
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(INCLUDE)

# AvgSbgrphFinder

AVGSBGRPHFINDER_INCLUDE_DEPS=include/deregnet/AvgSbgrphFinder.h \
                             include/deregnet/AvgSbgrphModel.h \
				             include/deregnet/utils.h \
				  		     include/deregnet/usinglemon.h

build/AvgSbgrphFinder.o : src/AvgSbgrphFinder.cpp $(AVGSBGRPHFINDER_INCLUDE_DEPS)
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(INCLUDE)

# AvgSbgrphModel

AVGSBGRPHMODEL_INCLUDE_DEPS=$(AVGSBGRPHFINDER_INCLUDE_DEPS)

build/AvgSbgrphModel.o : src/AvgSbgrphModel.cpp $(AVGSBGRPHMODEL_INCLUDE_DEPS)
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(INCLUDE)


clean : 
	rm -f $(objects)

destroy :
	rm -f $(objects)
	rm -f bin/avgdrgnt

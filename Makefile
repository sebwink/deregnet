SELF := $(abspath $(firstword $(MAKEFILE_LIST)))                                                                                                                          
SELFDIR := $(dir $(SELF))

GUROBI_MAKEFILE ?= $(SELFDIR)/upstream/libgrbfrc/gurobi.mak 

-include ${GUROBI_MAKEFILE}

CXX=c++
CXXFLAGS += -m64 -g -std=c++17 -Wall -Wextra -pedantic -fPIC #-Werror
CPPFLAGS += -Isrc -I${GUROBI_HOME}/include -I$(SELFDIR)/upstream/libgrbfrc/include
LIBGRBFRC ?= $(SELFDIR)/upstream/libgrbfrc/lib
LDPATHS += -L${GUROBI_HOME}/lib -lgurobi${GUROBI_VERSION_SUFFIX} -lgurobi_c++ -L$(LIBGRBFRC) -lgrbfrc
LDFLAGS += ${LDPATHS} -Wl,-rpath=$(LIBGRBFRC)
LDLIBS += -lpthread -lm

SRCDIR=src
SOURCES=$(wildcard $(SRCDIR)/*.cpp) 
HEADERS=$(wildcard $(SRCDIR)/*.hpp) 
OBJECTS=$(patsubst %.cpp, %.o, $(SOURCES))
DEPS=$(patsubst %.cpp, %.d, $(SOURCES))

BINDIR=src/bin
BIN_SOURCES=$(wildcard $(BINDIR)/*.cpp)
BINS=$(patsubst %.cpp, %, $(BIN_SOURCES)) 
BINDEPS=$(patsubst %.cpp, %.d, $(BIN_SOURCES))
BINOS=$(patsubst %.cpp, %.o, $(BIN_SOURCES)) 

clean:
	rm -f $(OBJECTS) $(DEPS) $(BINOS) $(BINDEPS)
	@cd $(SELFDIR)/upstream/libgrbfrc && make clean && rm -rf include

destroy: clean
	rm -rf doc/doxygen/html doc/doxygen/latex bin
	@cd $(SELFDIR)/upstream/libgrbfrc && make destroy

deregnet: install-python lgrbfrc all 
	@#


all : $(BINS)
	@#

$(BINS) : % : %.o 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)
	@mkdir -p bin && cp $@ bin/

$(BINS) : $(OBJECTS)

$(BINDEPS) : %.d : %.cpp
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -MM $< | sed 's,.*.o:,$*.o:,g' > $@

$(BINOS): %.o : %.cpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

lgrbfrc:
	@cd $(SELFDIR)/upstream/libgrbfrc && \
			make all && \
		    mkdir -p include/grbfrc && \
			cp src/*.hpp include/grbfrc

$(DEPS): %.d : %.cpp
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -MM $< | sed 's,.*.o:,$*.o:,g' > $@

$(OBJECTS): %.o : %.cpp  
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

docker: docker@$(GUROBI_VERSION)
	@#

docker@%:
	@make destroy
	@cd upstream/libgrbfrc && git submodule update --init --recursive && cd upstream/gurobi-docker && make gurobi@$*-local
	@cd upstream/libgrbfrc && make docker@$* 
	@docker build --no-cache --build-arg GUROBI_USER=$$(pwd) --build-arg GRBFRC_GUROBI_VERSION=$* --build-arg GRBFRC_IMAGE_TAG=local -t sebwink/deregnet:local .

docs:  
	@cd doc/doxygen && doxygen Doxyfile

docs-open: docs 
	@cd doc/doxygen/html && $(BROWSER) index.html

install-python:
	python3 -m pip install biomap-utils

-include $(DEPS)
-include $(BINDEPS)

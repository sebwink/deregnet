CXX = c++
CXXFLAGS = -std=c++11 -Wall -pedantic -fPIC #-ggdb3 -pg # -fopenmp # -DGRBFRC_OMP

include ../gurobi_version.mak

INCLUDE = -I${GUROBI_HOME}/include -Iinclude
LDFLAGS = -pg -L${GUROBI_HOME}/lib -lgurobi_c++ -lgurobi${GUROBI_VERSION_SUFFIX}

objects = FMILP.o CharnesCooper.o Dinkelbach.o Gloverizer.o YGGY.o Solver.o Algorithm.o ObjVarTransform.o

all : dirs libgrbfrc.a libgrbfrc.so 

doc : 
	doxygen doc/user/Doxyfile

devdoc :
	doxygen doc/dev/Doxyfile

dirs :
	mkdir -p lib

libgrbfrc.a : $(objects)
	ar cr $@ $(objects)
	mv $@ lib

FMILP.o : src/FMILP.cpp include/grbfrc/FMILP.h include/grbfrc/CharnesCooper.h include/grbfrc/Dinkelbach.h
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $<


%.o : src/%.cpp include/grbfrc/%.h
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $<

clean :
	rm -fr lib/libgrbfrc.a lib/libgrbfrc.so *.o

libgrbfrc.so : $(objects)
	$(CXX) $(objects) -o $@ -shared $(LDFLAGS) -Wl,-rpath=${GUROBI_HOME}/lib
	mv $@ lib
	rm *.o

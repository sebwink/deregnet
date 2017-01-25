CXX = c++
CXXFLAGS = -std=c++11 -Wall -pedantic -fPIC

GRB_PATH = /opt/gurobi/7.0.1

INCLUDE = -I${GRB_PATH}/include -Iinclude
LDFLAGS = -L${GRB_PATH}/lib -lgurobi_c++ -lgurobi70

objects = FMILP.o CharnesCooper.o Dinkelbach.o Gloverizer.o YGGY.o

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
	$(CXX) $(objects) -o $@ -shared $(LDFLAGS)
	mv $@ lib
	rm *.o

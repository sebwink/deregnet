EX= example1

GRB_HOME = /opt/gurobi/7.0.2
GRBFRC_HOME = ..

EXAMPLES = example1 example2 example3 example4 example5 example6 \
           example7 example8 example9

CXX = c++
CPPFLAGS = -std=c++11 -Wall -pedantic
INCLUDE_FLAGS = -I../include -I${GRB_HOME}/include 
LDFLAGS = -L${GRB_HOME}/lib -lgurobi_c++ -lgurobi70 -L../lib -lgrbfrc
RUNPATH = -Wl,-rpath=${GRBFRC_HOME}/lib

${EX} : ${EX}.o
	$(CXX) -o $@ $^ $(INCLUDE_FLAGS) $(LDFLAGS) $(RUNPATH)

${EX}.o : ${EX}.cpp
	$(CXX) $(CPPFLAGS) -c $< $(INCLUDE_FLAGS) 

superclean :
	rm -f *.o *.log $(EXAMPLES)

clean : 
	rm -f *.o *.log
